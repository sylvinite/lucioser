"""Creates sbatch jobs for MLST instances
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import glob
import gzip
import json
import os
import re
import shutil
import subprocess
import time

from datetime import datetime
from pathlib import Path

import yaml

from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.utils.referencer import Referencer


class Job_Creator:
    def __init__(self, config, log, sampleinfo={}, run_settings={}):
        self.config = config
        self.logger = log
        self.batchfile = "/tmp/batchfile.sbatch"

        self.filelist = list()
        if isinstance(run_settings.get("input"), list):
            self.filelist = run_settings.get("input")
            run_settings["input"] = "/tmp/"

        self.run_settings = run_settings
        self.indir = os.path.abspath(run_settings.get("input", "/tmp/"))
        self.trimmed = run_settings.get("trimmed", True)
        self.qc_only = run_settings.get("qc_only", False)
        self.careful = run_settings.get("careful", True)
        self.pool = run_settings.get("pool", [])
        self.finishdir = run_settings.get("finishdir", "")

        self.sampleinfo = sampleinfo
        self.sample = None
        if isinstance(self.sampleinfo, list) and len(self.sampleinfo) > 1:
            self.name = self.sampleinfo[0].get("CG_ID_project")
            self.sample = self.sampleinfo[0]
            for entry in self.sampleinfo:
                if entry.get("CG_ID_sample") == self.name:
                    raise Exception(
                        "Mixed projects in samples_info file. Do not know how to proceed"
                    )
        else:
            if isinstance(self.sampleinfo, list):
                self.sampleinfo = self.sampleinfo[0]
            self.name = self.sampleinfo.get("CG_ID_sample")
            self.sample = self.sampleinfo

        # If timestamp is provided. Use it as analysis time. Else use current time
        if run_settings.get("timestamp") is not None:
            self.now = run_settings.get("timestamp")
            temp = run_settings.get("timestamp").replace("_", ".").split(".")
            self.dt = datetime(
                int(temp[0]),
                int(temp[1]),
                int(temp[2]),
                int(temp[3]),
                int(temp[4]),
                int(temp[5]),
            )
        else:
            self.dt = datetime.now()
            self.now = time.strftime(
                "{}.{}.{}_{}.{}.{}".format(
                    self.dt.year,
                    self.dt.month,
                    self.dt.day,
                    self.dt.hour,
                    self.dt.minute,
                    self.dt.second,
                )
            )

        if run_settings.get("finishdir") is None:
            self.finishdir = "{}/{}_{}".format(config["folders"]["results"], self.name, self.now)
        self.db_pusher = DB_Manipulator(config, log)
        self.concat_files = dict()
        self.ref_resolver = Referencer(config, log)

    def get_sbatch(self):
        """ Returns sbatchfile, slightly superflous"""
        return self.batchfile

    def get_headerargs(self):
        headerline = "-A {} -p {} -n {} -t {} -J {}_{} --qos {} --output {}/slurm_{}.log".format(
            self.config["slurm_header"]["project"],
            self.config["slurm_header"]["type"],
            self.config["slurm_header"]["threads"],
            self.config["slurm_header"]["time"],
            self.config["slurm_header"]["job_prefix"],
            self.name,
            self.config["slurm_header"]["qos"],
            self.finishdir,
            self.name,
        )
        return headerline

    def verify_fastq(self):
        """ Uses arg indir to return a dict of PE fastq tuples fulfilling naming convention """
        verified_files = list()
        files = os.listdir(self.indir)
        if files == []:
            raise Exception("Directory {} lacks fastq files.".format(self.indir))
        for file in files:
            file_match = re.match(self.config["regex"]["file_pattern"], file)
            if file_match:
                # Check that symlinks resolve
                path = "{}/{}".format(self.indir, file)
                if os.path.islink(path):
                    if not os.path.exists(os.readlink(path)):
                        raise Exception(
                            "Some fastq files are unresolved symlinks in directory {}.".format(
                                self.indir
                            )
                        )

                # Make sure both mates exist
                if (
                    file_match[1] == "1"
                    or file_match[1] == "2"
                    or file_match[1] == "forward"
                    or file_match[1] == "reverse"
                ):
                    if file_match[1] == "forward" or file_match[1] == "reverse":
                        pairno = "forward"
                        if "forward" in file_match[1]:
                            pairno = "reverse"
                        pairname = file_match[0].replace(file_match[1], pairno)
                    else:
                        pairno = 2 - 1 % int(file_match[1])  # 1->2, 2->1
                        # Construct mate name
                        pairname = "{}{}{}".format(
                            file_match.string[: file_match.end(1) - 1],
                            pairno,
                            file_match.string[file_match.end(1) : file_match.end()],
                        )
                    if pairname in files:
                        files.pop(files.index(pairname))
                        verified_files.append(file_match[0])
                        verified_files.append(pairname)
                else:
                    raise Exception(
                        "Some fastq files have no mate in directory {}.".format(self.indir)
                    )
        if verified_files == []:
            raise Exception(
                "No files in directory {} match file_pattern '{}'.".format(
                    self.indir, self.config["regex"]["file_pattern"]
                )
            )

        # Warn about file sizes
        for vfile in verified_files:
            try:
                bsize = os.stat("{}/{}".format(self.indir, vfile)).st_size
                bsize = bsize >> 20
                if bsize > 1000:
                    self.logger.warning("Input fastq {} exceeds 1000MB".format(vfile))
            except Exception as e:
                self.logger.warning(
                    "Unable to verify size of input file {}/{}".format(self.indir, vfile)
                )

        # Warn about invalid fastq files
        for vfile in verified_files:
            f = gzip.open("{}/{}".format(self.indir, vfile), "r")
            lines = f.read().splitlines()
            if len(lines) < 2 or not "+" in str(lines[-2]):
                self.logger.warning("Input fastq {} does not seem to end properly".format(vfile))
        return sorted(verified_files)

    def create_assemblysection(self):
        batchfile = open(self.batchfile, "a+")
        # memory is actually 128 per node regardless of cores.
        batchfile.write("# Spades assembly\n")
        if self.trimmed:
            trimline = "-s {}".format(self.concat_files["i"])
        else:
            trimline = ""
        if self.careful:
            careline = "--careful"
        else:
            careline = ""

        batchfile.write(
            "spades.py --threads {} {} --memory {} -o {}/assembly -1 {} -2 {} {}\n".format(
                self.config["slurm_header"]["threads"],
                careline,
                8 * int(self.config["slurm_header"]["threads"]),
                self.finishdir,
                self.concat_files["f"],
                self.concat_files["r"],
                trimline,
            )
        )

        batchfile.write(
            "mv {0}/assembly/contigs.fasta {0}/assembly/{1}_contigs.fasta\n".format(
                self.finishdir, self.name
            )
        )
        batchfile.write(
            "sed -n '/NODE_1000_/q;p' {0}/assembly/{1}_contigs.fasta > {0}/assembly/{1}_trimmed_contigs.fasta\n".format(
                self.finishdir, self.name
            )
        )
        # batchfile.write("##Input cleanup\n")
        # batchfile.write("rm -r {}/trimmed\n".format(self.finishdir))
        batchfile.write("\n\n")
        batchfile.close()

    def blast_subset(self, name, search_string):
        # Create run
        file_list = glob.glob(search_string)
        batchfile = open(self.batchfile, "a+")
        batchfile.write("mkdir {}/blast_search/{}\n".format(self.finishdir, name))
        blast_format = '"7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send length"'

        if len(file_list) > 1:
            for ref in file_list:
                if re.search(r"(\w+(?:\-\w+)*)\.\w+", os.path.basename(ref)) is None:
                    self.logger.error(
                        "File {} does not match typical format. Consider deleting and redownloading"
                    )
                else:
                    ref_nosuf = re.search(r"(\w+(?:\-\w+)*)\.\w+", os.path.basename(ref)).group(1)
                batchfile.write(
                    "# BLAST {} search for {}, {}\n".format(
                        name, self.sample.get("organism"), ref_nosuf
                    )
                )
                if name == "mlst":
                    batchfile.write(
                        "blastn -db {}/{}  -query {}/assembly/{}_contigs.fasta -out {}/blast_search/{}/loci_query_{}.txt -task megablast -num_threads {} -outfmt {}\n".format(
                            os.path.dirname(ref),
                            ref_nosuf,
                            self.finishdir,
                            self.name,
                            self.finishdir,
                            name,
                            ref_nosuf,
                            self.config["slurm_header"]["threads"],
                            blast_format,
                        )
                    )
                else:
                    batchfile.write(
                        "blastn -db {}/{}  -query {}/assembly/{}_contigs.fasta -out {}/blast_search/{}/{}.txt -task megablast -num_threads {} -outfmt {}\n".format(
                            os.path.dirname(ref),
                            ref_nosuf,
                            self.finishdir,
                            self.name,
                            self.finishdir,
                            name,
                            ref_nosuf,
                            self.config["slurm_header"]["threads"],
                            blast_format,
                        )
                    )
        elif len(file_list) == 1:
            ref_nosuf = re.search(r"(\w+(?:\-\w+)*)\.\w+", os.path.basename(file_list[0])).group(1)
            batchfile.write(
                "## BLAST {} search in {}\n".format(
                    name, self.sample.get("organism").replace("_", " ").capitalize()
                )
            )
            batchfile.write(
                "blastn -db {}/{}  -query {}/assembly/{}_contigs.fasta -out {}/blast_search/{}/{}.txt -task megablast -num_threads {} -outfmt {}\n".format(
                    os.path.dirname(search_string),
                    ref_nosuf,
                    self.finishdir,
                    self.name,
                    self.finishdir,
                    name,
                    ref_nosuf,
                    self.config["slurm_header"]["threads"],
                    blast_format,
                )
            )
        batchfile.write("\n")
        batchfile.close()

    def create_variantsection(self):
        """ Creates a job for variant calling based on local alignment """
        ref = "{}/{}.fasta".format(self.config["folders"]["genomes"], self.sample.get("reference"))
        localdir = "{}/alignment".format(self.finishdir)
        outbase = "{}/{}_{}".format(localdir, self.name, self.sample.get("reference"))

        # Create run
        batchfile = open(self.batchfile, "a+")
        batchfile.write("# Variant calling based on local alignment\n")
        batchfile.write("mkdir {}\n".format(localdir))

        batchfile.write("## Alignment & Deduplication\n")
        batchfile.write(
            "bwa mem -M -t {} {} {} {} > {}.sam\n".format(
                self.config["slurm_header"]["threads"],
                ref,
                self.concat_files["f"],
                self.concat_files["r"],
                outbase,
            )
        )
        batchfile.write(
            "samtools view --threads {} -b -o {}.bam -T {} {}.sam\n".format(
                self.config["slurm_header"]["threads"], outbase, ref, outbase
            )
        )
        batchfile.write(
            "samtools sort --threads {} -o {}.bam_sort {}.bam\n".format(
                self.config["slurm_header"]["threads"], outbase, outbase
            )
        )
        batchfile.write(
            "picard MarkDuplicates I={}.bam_sort O={}.bam_sort_rmdup M={}.stats.dup REMOVE_DUPLICATES=true\n".format(
                outbase, outbase, outbase
            )
        )
        batchfile.write("samtools index {}.bam_sort_rmdup\n".format(outbase))
        batchfile.write(
            "samtools idxstats {}.bam_sort_rmdup &> {}.stats.ref\n".format(outbase, outbase)
        )
        # Removal of temp aligment files
        batchfile.write("rm {}.bam {}.sam\n".format(outbase, outbase))

        batchfile.write("## Primary stats generation\n")
        # Insert stats, dedupped
        batchfile.write(
            "picard CollectInsertSizeMetrics I={}.bam_sort_rmdup O={}.stats.ins H={}.hist.ins\n".format(
                outbase, outbase, outbase
            )
        )
        # Coverage
        batchfile.write(
            "samtools stats --coverage 1,10000,1 {}.bam_sort_rmdup |grep ^COV | cut -f 2- &> {}.stats.cov\n".format(
                outbase, outbase
            )
        )
        # Mapped rate, no dedup,dedup in MWGS (trimming has no effect)!
        batchfile.write("samtools flagstat {}.bam_sort &> {}.stats.map\n".format(outbase, outbase))
        # Total reads, no dedup,dedup in MWGS (trimming has no effect)!
        batchfile.write("samtools view -c {}.bam_sort &> {}.stats.raw\n".format(outbase, outbase))

        batchfile.write("\n\n")
        batchfile.close()

    def create_preprocsection(self):
        """Concatinates data, possibly trims it, then makes the unstranded reads usable"""
        forward = list()
        reverse = list()
        for root, dirs, files in os.walk(self.config["folders"]["adapters"]):
            if not "NexteraPE-PE.fa" in files:
                self.logger.error(
                    "Adapters folder at {} does not contain NexteraPE-PE.fa. Review paths.yml"
                )
            else:
                break
        trimdir = "{}/trimmed".format(self.finishdir)
        files = self.verify_fastq()
        batchfile = open(self.batchfile, "a+")
        batchfile.write("#Trimmomatic section\n")
        batchfile.write("mkdir {}\n".format(trimdir))

        batchfile.write("##Pre-concatination\n")
        for file in files:
            fullfile = "{}/{}".format(self.indir, file)
            # Even indexes = Forward
            if not files.index(file) % 2:
                forward.append(fullfile)
            elif files.index(file) % 2:
                reverse.append(fullfile)
        outfile = files[0].split("_")[0]

        self.concat_files["f"] = "{}/trimmed/{}_forward_reads.fastq.gz".format(
            self.finishdir, self.name
        )
        self.concat_files["r"] = "{}/trimmed/{}_reverse_reads.fastq.gz".format(
            self.finishdir, self.name
        )
        batchfile.write("cat {} > {}\n".format(" ".join(forward), self.concat_files.get("f")))
        batchfile.write("cat {} > {}\n".format(" ".join(reverse), self.concat_files.get("r")))

        if self.trimmed:
            fp = "{}/{}_trim_front_pair.fastq.gz".format(trimdir, outfile)
            fu = "{}/{}_trim_front_unpair.fastq.gz".format(trimdir, outfile)
            rp = "{}/{}_trim_rev_pair.fastq.gz".format(trimdir, outfile)
            ru = "{}/{}_trim_rev_unpair.fastq.gz".format(trimdir, outfile)
            batchfile.write("##Trimming section\n")
            batchfile.write(
                "trimmomatic PE -threads {} -phred33 {} {} {} {} {} {}\
      ILLUMINACLIP:{}/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\n".format(
                    self.config["slurm_header"]["threads"],
                    self.concat_files.get("f"),
                    self.concat_files.get("r"),
                    fp,
                    fu,
                    rp,
                    ru,
                    self.config["folders"]["adapters"],
                )
            )

            batchfile.write("## Interlaced trimmed files\n")
            self.concat_files["f"] = fp
            self.concat_files["r"] = rp
            self.concat_files["i"] = "{}/{}_trim_unpair.fastq.gz".format(trimdir, outfile)

            batchfile.write("cat {} >> {}\n".format(" ".join([fu, ru]), self.concat_files.get("i")))
        batchfile.write("\n")
        batchfile.close()

    def create_assemblystats_section(self):
        batchfile = open(self.batchfile, "a+")
        batchfile.write("# QUAST QC metrics\n")
        batchfile.write("mkdir {}/assembly/quast\n".format(self.finishdir))
        batchfile.write(
            "quast.py {}/assembly/{}_contigs.fasta -o {}/assembly/quast\n".format(
                self.finishdir, self.name, self.finishdir
            )
        )
        batchfile.write(
            "mv {}/assembly/quast/report.tsv {}/assembly/quast/{}_report.tsv\n\n".format(
                self.finishdir, self.finishdir, self.name
            )
        )
        batchfile.close()

    def create_snpsection(self):
        snplist = self.filelist.copy()
        batchfile = open(self.batchfile, "a+")
        name = ""

        # VCFTools filters:
        vcffilter = "--minQ 30 --thin 50 --minDP 3 --min-meanDP 20"
        # BCFTools filters:
        bcffilter = "GL[0]<-500 & GL[1]=0 & QR/RO>30 & QA/AO>30 & QUAL>5000 & ODDS>1100 & GQ>140 & DP>100 & MQM>59 & SAP<15 & PAIRED>0.9 & EPP>3"

        for item in snplist:
            if item.count("/") >= 2:
                name = item.split("/")[-2]
            if "_" in name:
                name = name.split("_")[0]
            batchfile.write("# Basecalling for sample {}\n".format(name))
            ref = "{}/{}.fasta".format(
                self.config["folders"]["genomes"], self.sample.get("reference")
            )
            outbase = "{}/{}_{}".format(item, name, self.sample.get("reference"))
            batchfile.write(
                "samtools view -h -q 1 -F 4 -F 256 {}.bam_sort_rmdup | grep -v XA:Z | grep -v SA:Z| samtools view -b - > {}/{}.unique\n".format(
                    outbase, self.finishdir, name
                )
            )
            batchfile.write(
                "freebayes -= --pvar 0.7 -j -J --standard-filters -C 6 --min-coverage 30 --ploidy 1 -f {} -b {}/{}.unique -v {}/{}.vcf\n".format(
                    ref, self.finishdir, name, self.finishdir, name
                )
            )
            batchfile.write(
                "bcftools view {}/{}.vcf -o {}/{}.bcf.gz -O b --exclude-uncalled --types snps\n".format(
                    self.finishdir, name, self.finishdir, name
                )
            )
            batchfile.write("bcftools index {}/{}.bcf.gz\n".format(self.finishdir, name))
            batchfile.write("\n")

            batchfile.write(
                "vcftools --bcf {}/{}.bcf.gz {} --remove-filtered-all --recode-INFO-all --recode-bcf --out {}/{}\n".format(
                    self.finishdir, name, vcffilter, self.finishdir, name
                )
            )
            batchfile.write(
                'bcftools view {}/{}.recode.bcf -i "{}" -o {}/{}.recode.bcf.gz -O b --exclude-uncalled --types snps\n'.format(
                    self.finishdir, name, bcffilter, self.finishdir, name
                )
            )
            batchfile.write("bcftools index {}/{}.recode.bcf.gz\n\n".format(self.finishdir, name))

        batchfile.write("# SNP pair-wise distance\n")
        batchfile.write("touch {}/stats.out\n".format(self.finishdir))
        while len(snplist) > 1:
            nameOne = ""
            nameTwo = ""
            top = snplist.pop(0)
            if top.count("/") >= 2:
                nameOne = top.split("/")[-2]
            if "_" in nameOne:
                nameOne = nameOne.split("_")[0]
            for entry in snplist:
                if entry.count("/") >= 2:
                    nameTwo = entry.split("/")[-2]
                if "_" in nameTwo:
                    nameTwo = nameTwo.split("_")[0]

                pair = "{}_{}".format(nameOne, nameTwo)
                batchfile.write(
                    "bcftools isec {}/{}.recode.bcf.gz {}/{}.recode.bcf.gz -n=1 -c all -p {}/tmp -O b\n".format(
                        self.finishdir, nameOne, self.finishdir, nameTwo, self.finishdir
                    )
                )
                batchfile.write(
                    "bcftools merge -O b -o {}/{}.bcf.gz --force-samples {}/tmp/0000.bcf {}/tmp/0001.bcf\n".format(
                        self.finishdir, pair, self.finishdir, self.finishdir
                    )
                )
                batchfile.write("bcftools index {}/{}.bcf.gz\n".format(self.finishdir, pair))

                batchfile.write(
                    "echo {} $( bcftools stats {}/{}.bcf.gz |grep SNPs: | cut -d $'\\t' -f4 ) >> {}/stats.out\n".format(
                        pair, self.finishdir, pair, self.finishdir
                    )
                )
                batchfile.write("\n")
        batchfile.close()

    def create_collection(self):
        """Creates collection entry in database"""
        if self.db_pusher.exists("Collections", {"ID_collection": self.name}):
            self.db_pusher.purge_rec(name=self.name, type="Collections")
            for sample in self.pool:
                self.db_pusher.add_rec(
                    {"ID_collection": self.name, "CG_ID_sample": sample}, "Collections"
                )

        addedprojs = list()
        for sample in self.pool:
            proj = re.search(r"(\w+)A(?:\w+)", sample).group(1)
            if proj not in addedprojs:
                self.create_project(proj)
                addedprojs.append(proj)

    def create_project(self, name):
        """Creates project in database"""
        proj_col = dict()
        proj_col["CG_ID_project"] = name
        proj_col["Customer_ID_project"] = self.sample.get("Customer_ID_project")
        proj_col["Customer_ID"] = self.sample.get("Customer_ID")
        self.db_pusher.add_rec(proj_col, "Projects")
        self.db_pusher.upd_rec({"CG_ID_project": name}, "Projects", proj_col)

    def create_sample(self, name):
        """Creates sample in database"""
        try:
            sample_col = self.db_pusher.get_columns("Samples")
            sample_col["CG_ID_sample"] = self.sample.get("CG_ID_sample")
            sample_col["CG_ID_project"] = self.sample.get("CG_ID_project")
            sample_col["Customer_ID_sample"] = self.sample.get("Customer_ID_sample")
            sample_col["reference_genome"] = self.sample.get("reference")
            sample_col["date_analysis"] = self.dt
            sample_col["organism"] = self.sample.get("organism")
            sample_col["application_tag"] = self.sample.get("application_tag")
            sample_col["priority"] = self.sample.get("priority")
            sample_col["date_arrival"] = datetime.strptime(
                self.sample.get("date_arrival"), "%Y-%m-%d %H:%M:%S"
            )
            sample_col["date_sequencing"] = datetime.strptime(
                self.sample.get("date_sequencing"), "%Y-%m-%d %H:%M:%S"
            )
            sample_col["date_libprep"] = datetime.strptime(
                self.sample.get("date_libprep"), "%Y-%m-%d %H:%M:%S"
            )
            sample_col["method_libprep"] = self.sample.get("method_libprep")
            sample_col["method_sequencing"] = self.sample.get("method_sequencing")
            # self.db_pusher.purge_rec(sample_col['CG_ID_sample'], 'sample')
            self.db_pusher.add_rec(sample_col, "Samples")
        except Exception as e:
            self.logger.error("Unable to add sample {} to database".format(self.name))

    def project_job(self, single_sample=False):
        if "dry" in self.config and self.config["dry"] == True:
            dry = True
        else:
            dry = False
        jobarray = list()
        if not os.path.exists(self.finishdir):
            os.makedirs(self.finishdir)
        # Loads project level info.
        try:
            if single_sample:
                self.create_project(self.sample.get("CG_ID_project"))
            elif self.pool:
                self.create_collection()
            else:
                self.create_project(self.name)
        except Exception as e:
            self.logger.error(
                "LIMS interaction failed. Unable to read/write project {}".format(self.name)
            )
        # Writes the job creation sbatch
        if single_sample:
            try:
                self.sample_job()
                headerargs = self.get_headerargs()
                outfile = self.get_sbatch()
                bash_cmd = "sbatch {} {}".format(headerargs, outfile)
                if not dry and outfile != "":
                    samproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
                    output, error = samproc.communicate()
                    jobno = re.search(r"(\d+)", str(output)).group(0)
                    jobarray.append(jobno)
                else:
                    self.logger.info("Suppressed command: {}".format(bash_cmd))
            except Exception as e:
                self.logger.error("Unable to analyze single sample {}".format(self.name))
        else:
            for ldir in glob.glob("{}/*/".format(self.indir)):
                ldir = os.path.basename(os.path.normpath(ldir))
                try:
                    sample_in = "{}/{}".format(self.indir, ldir)
                    sample_out = "{}/{}".format(self.finishdir, ldir)
                    linkedjson = None
                    local_sampleinfo = [p for p in self.sampleinfo if p["CG_ID_sample"] == ldir]
                    if local_sampleinfo == []:
                        raise Exception("Sample {} has no counterpart in json file".format(ldir))
                    else:
                        local_sampleinfo = local_sampleinfo[0]
                    sample_settings = dict(self.run_settings)
                    sample_settings["input"] = sample_in
                    sample_settings["finishdir"] = sample_out
                    sample_settings["timestamp"] = self.now
                    sample_instance = Job_Creator(
                        config=self.config,
                        log=self.logger,
                        sampleinfo=local_sampleinfo,
                        run_settings=sample_settings,
                    )
                    sample_instance.sample_job()
                    headerargs = sample_instance.get_headerargs()
                    outfile = ""
                    if os.path.isfile(sample_instance.get_sbatch()):
                        outfile = sample_instance.get_sbatch()
                    bash_cmd = "sbatch {} {}".format(headerargs, outfile)
                    if not dry and outfile != "":
                        projproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
                        output, error = projproc.communicate()
                        jobno = re.search(r"(\d+)", str(output)).group(0)
                        jobarray.append(jobno)
                    else:
                        self.logger.info("Suppressed command: {}".format(bash_cmd))
                except Exception as e:
                    pass
        if not dry:
            self.finish_job(jobarray, single_sample)

    def finish_job(self, joblist, single_sample=False):
        """ Uploads data and sends an email once all analysis jobs are complete. """
        report = "default"
        if self.qc_only:
            report = "qc"
        custom_conf = ""
        if "config_path" in self.config:
            custom_conf = "--config {}".format(self.config["config_path"])

        process = subprocess.Popen("id -un".split(), stdout=subprocess.PIPE)
        user, error = process.communicate()
        user = str(user).replace(".", " ").title()
        # if not os.path.exists(self.finishdir):
        #  os.makedirs(self.finishdir)

        startfile = "{}/run_started.out".format(self.finishdir)
        configfile = "{}/config.log".format(self.finishdir)
        mailfile = "{}/mailjob.sh".format(self.finishdir)
        samplefile = "{}/sampleinfo.json".format(self.finishdir)
        with open(samplefile, "w+") as outfile:
            json.dump(self.sampleinfo, outfile)

        sb = open(startfile, "w+")
        cb = open(configfile, "w+")
        mb = open(mailfile, "w+")

        sb.write("#!/usr/bin/env bash\n")
        sb.close()
        configout = self.config.copy()
        if "genologics" in configout:
            del configout["genologics"]
        cb.write("ANALYSIS STARTED BY: {}\n".format(user))
        cb.write(json.dumps(configout, indent=2, separators=(",", ":")))
        cb.close()
        mb.write("#!/usr/bin/env bash\n\n")
        mb.write("#Uploading of results to database and production of report\n")
        if "MICROSALT_CONFIG" in os.environ:
            mb.write("export MICROSALT_CONFIG={}\n".format(os.environ["MICROSALT_CONFIG"]))
        mb.write("source activate $CONDA_DEFAULT_ENV\n")

        ct = ""
        if "custom_target" in self.run_settings and self.run_settings["custom_target"] != "":
            ct = "--custom_target {}".format(self.run_settings["custom_target"])

        mb.write(
            "microSALT utils finish {0}/sampleinfo.json --input {0} --email {1} --report {2} {3} {4}\n".format(
                self.finishdir,
                self.config["regex"]["mail_recipient"],
                report,
                custom_conf,
                ct,
            )
        )
        mb.write("touch {}/run_complete.out".format(self.finishdir))
        mb.close()

        massagedJobs = list()
        final = ":".join(joblist)
        # Create subtracker if more than 50 samples
        maxlen = 50
        if len(joblist) > maxlen:
            i = 1
            while i <= len(joblist):
                if i + maxlen < len(joblist):
                    massagedJobs.append(":".join(joblist[i - 1 : i + maxlen - 1]))
                else:
                    massagedJobs.append(":".join(joblist[i - 1 : -1]))
                i += maxlen
            for entry in massagedJobs:
                if massagedJobs.index(entry) < len(massagedJobs) - 1:
                    head = "-A {} -p core -n 1 -t 00:00:10 -J {}_{}_SUBTRACKER --qos {} --dependency=afterany:{}".format(
                        self.config["slurm_header"]["project"],
                        self.config["slurm_header"]["job_prefix"],
                        self.name,
                        self.config["slurm_header"]["qos"],
                        entry,
                    )
                    bash_cmd = "sbatch {} {}".format(head, startfile)
                    mailproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
                    output, error = mailproc.communicate()
                    jobno = re.search(r"(\d+)", str(output)).group(0)
                    massagedJobs[massagedJobs.index(entry) + 1] += ":{}".format(jobno)
                else:
                    final = entry
                    break

        head = "-A {} -p core -n 1 -t 6:00:00 -J {}_{}_MAILJOB --qos {} --open-mode append --dependency=afterany:{} --output {}".format(
            self.config["slurm_header"]["project"],
            self.config["slurm_header"]["job_prefix"],
            self.name,
            self.config["slurm_header"]["qos"],
            final,
            self.config["folders"]["log_file"],
            self.config["regex"]["mail_recipient"],
        )
        bash_cmd = "sbatch {} {}".format(head, mailfile)
        mailproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
        output, error = mailproc.communicate()

        try:
            jobno = str(re.search(r"(\d+)", str(output)).group(0))
            joblist.append(jobno)
        except Exception as e:
            self.logger.info("Unable to grab SLURMID for {0}".format(self.name))

        try:
            #Generates file with all slurm ids
            slurmname = "{}_slurm_ids.yaml".format(self.name)
            slurmreport_storedir = Path(self.config["folders"]["reports"],
                "trailblazer", slurmname)
            slurmreport_workdir = Path(self.finishdir, slurmname)
            yaml.safe_dump(
                data={"jobs": [str(job) for job in joblist]},
                      stream=open(slurmreport_workdir, "w"))
            shutil.copyfile(slurmreport_workdir, slurmreport_storedir)
            self.logger.info(
                "Saved Trailblazer slurm report file to %s and %s",
                slurmreport_storedir,
                slurmreport_workdir,
            )
        except Exception as e:
            self.logger.info("Unable to generate Trailblazer slurm report file")

    def sample_job(self):
        """ Writes necessary sbatch job for each individual sample """
        try:
            if not os.path.exists(self.finishdir):
                os.makedirs(self.finishdir)
            try:
                # This is one job
                self.batchfile = "{}/runfile.sbatch".format(self.finishdir)
                batchfile = open(self.batchfile, "w+")
                batchfile.write("#!/bin/sh\n\n")
                batchfile.write("mkdir -p {}\n".format(self.finishdir))
                batchfile.close()
                self.create_preprocsection()
                self.create_variantsection()
                if not self.qc_only:
                    self.create_assemblysection()
                    self.create_assemblystats_section()
                    self.create_blast_search()
                batchfile = open(self.batchfile, "a+")
                batchfile.close()

                self.logger.info(
                    "Created runfile for sample {} in folder {}".format(self.name, self.finishdir)
                )
            except Exception as e:
                raise
            try:
                self.create_sample(self.name)
            except Exception as e:
                self.logger.error("Unable to access LIMS info for sample {}".format(self.name))
        except Exception as e:
            self.logger.error(
                "Unable to create job for sample {}\nSource: {}".format(self.name, str(e))
            )
            shutil.rmtree(self.finishdir, ignore_errors=True)
            raise

    def create_blast_search(self):
        reforganism = self.ref_resolver.organism2reference(self.sample.get("organism"))
        self.batchfile = "{}/runfile.sbatch".format(self.finishdir)
        batchfile = open(self.batchfile, "a+")
        batchfile.write("mkdir -p {}/blast_search\n".format(self.finishdir))
        batchfile.close()
        self.blast_subset(
            "mlst",
            "{}/{}/*.tfa".format(self.config["folders"]["references"], reforganism),
        )
        
        self.blast_subset(
            "resistance", "{}/*.fsa".format(self.config["folders"]["resistances"])
        )
        if self.run_settings["custom_target"] != "":
            self.blast_subset("custom", "{}/*.fsa".format(self.run_settings["custom_target"]))


    def snp_job(self):
        """ Writes a SNP calling job for a set of samples """
        if not os.path.exists(self.finishdir):
            os.makedirs(self.finishdir)

        self.batchfile = "{}/runfile.sbatch".format(self.finishdir)
        batchfile = open(self.batchfile, "w+")
        batchfile.write("#!/usr/bin/env bash\n")
        batchfile.write("mkdir -p {}\n\n".format(self.finishdir))
        batchfile.close()

        self.create_snpsection()
        batchfile = open(self.batchfile, "a+")
        batchfile.close()

        headerline = (
            "-A {} -p {} -n 1 -t 24:00:00 -J {}_{} --qos {} --output {}/slurm_{}.log".format(
                self.config["slurm_header"]["project"],
                self.config["slurm_header"]["type"],
                self.config["slurm_header"]["job_prefix"],
                self.name,
                self.config["slurm_header"]["qos"],
                self.finishdir,
                self.name,
            )
        )
        outfile = self.get_sbatch()
        bash_cmd = "sbatch {} {}".format(headerline, outfile)
        samproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
        output, error = samproc.communicate()
