"""Creates sbatch jobs for MLST instances
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import glob
import os
import re
import shutil
import subprocess
import time

from datetime import datetime
from microSALT.store.lims_fetcher import LIMS_Fetcher
from microSALT.store.db_manipulator import DB_Manipulator

class Job_Creator():

  def __init__(self, input, config, log, finishdir="", timestamp=""):
    self.config = config
    self.logger = log
    self.batchfile = ""
    self.filelist = list()
    self.indir = ""

    if isinstance(input, str):
      self.indir = os.path.abspath(input)
      self.name = os.path.basename(os.path.normpath(self.indir))
    elif type(input) == list:
      self.filelist = input
      self.name = "SNP"

    self.now = timestamp
    if timestamp != "":
      self.now = timestamp
      temp = timestamp.replace('_','.').split('.')
      self.dt = datetime(int(temp[0]),int(temp[1]),int(temp[2]),int(temp[3]),int(temp[4]),int(temp[5]))
    else:
      self.dt = datetime.now() 
      self.now = time.strftime("{}.{}.{}_{}.{}.{}".\
      format(self.dt.year, self.dt.month, self.dt.day, self.dt.hour, self.dt.minute, self.dt.second))

    self.outdir = "/scratch/$SLURM_JOB_ID/workdir/{}_{}".format(self.name, self.now)
    self.finishdir = finishdir
    if self.finishdir == "":
      self.finishdir="{}/{}_{}".format(config["folders"]["results"], self.name, self.now)
    
    self.db_pusher=DB_Manipulator(config, log)
    self.trimmed_files = dict()
    self.concat_files = dict()
    self.organism = ""
    self.lims_fetcher = LIMS_Fetcher(config, log)

  def get_sbatch(self):
    """ Returns sbatchfile, slightly superflous"""
    return self.batchfile

  def get_headerargs(self):
    headerline = "-A {} -p {} -n {} -t {} -J {}_{} --qos {} --output {}/slurm_{}.log".format(self.config["slurm_header"]["project"],\
                 self.config["slurm_header"]["type"], self.config["slurm_header"]["threads"],self.config["slurm_header"]["time"],\
                 self.config["slurm_header"]["job_prefix"], self.name,self.config["slurm_header"]["qos"], self.finishdir, self.name)
    return headerline

  def verify_fastq(self):
    """ Uses arg indir to return a dict of PE fastq tuples fulfilling naming convention """
    files = os.listdir(self.indir)
    if files == []:
      raise Exception("Directory {} lacks fastq files.".format(self.indir))
    verified_files = dict()
    for file in files:
      file_match = re.match( self.config['regex']['file_pattern'], file)
      if file_match:
        #Check that symlinks resolve
        path = '{}/{}'.format(self.indir, file)
        if os.path.islink(path):
          if not os.path.exists(os.readlink(path)):
            raise Exception("Some fastq files are unresolved symlinks in directory {}.".format(self.indir))

        #Make sure both mates exist
        if file_match[1] == '1':
          pairno = '2'
          #Construct mate name
          pairname = "{}{}{}".format(file_match.string[:file_match.end(1)-1] , pairno, \
                      file_match.string[file_match.end(1):file_match.end()])
          if pairname in files:
            if not self.name in verified_files.keys():
              verified_files[self.name] = list()
            files.pop( files.index(pairname) )
            verified_files[self.name].append(file_match[0])
            verified_files[self.name].append(pairname)
        elif file_match[1] == '2':
          pass
        else:
          raise Exception("Some fastq files have no mate in directory {}.".format(self.indir))
    if verified_files == []:
      raise Exception("No files in directory {} match file_pattern '{}'.".format(self.indir, self.config['regex']['file_pattern']))
    return verified_files
 
  def interlace_files(self):
    """Interlaces all trimmed files"""
    fplist = list()
    kplist = list()
    ilist = list()
    batchfile = open(self.batchfile, "a+")
    batchfile.write("# Interlaced trimmed files\n")
 
    for name, v in self.trimmed_files.items():
      fplist.append( v['fp'] )
      kplist.append( v['rp'] )
      ilist.append( v['fu'] )
      ilist.append( v['ru'] )

    if len(kplist) != len(fplist) or len(ilist)/2 != len(kplist):
      raise Exception("Uneven distribution of trimmed files. Invalid trimming step {}".format(name))
    self.concat_files['f'] = "{}/trimmed/{}{}".format(self.outdir,self.name, "_trim_front_pair.fq")
    self.concat_files['r'] = "{}/trimmed/{}{}".format(self.outdir,self.name, "_trim_rev_pair.fq")
    self.concat_files['i'] = "{}/trimmed/{}{}".format(self.outdir,self.name, "_trim_unpaired.fq")
    for k, v in self.concat_files.items():
      batchfile.write("touch {}\n".format(v))
    
    batchfile.write("cat {} >> {}\n".format(' '.join(fplist), self.concat_files['f']))
    batchfile.write("cat {} >> {}\n".format(' '.join(kplist), self.concat_files['r']))
    batchfile.write("cat {} >> {}\n".format(' '.join(ilist), self.concat_files['i']))
    batchfile.write("rm {} {} {}\n".format(' '.join(fplist), ' '.join(kplist), ' '.join(ilist)))
    batchfile.write("\n")
    batchfile.close()    

  def create_assemblysection(self):
    batchfile = open(self.batchfile, "a+")
    #memory is actually 128 per node regardless of cores.
    batchfile.write("# Spades assembly\n")
    batchfile.write("cat {} > {}/trimmed/forward_reads.fasta.gz").format(self.concat_files['f'], self.outdir))
    batchfile.write("cat {} > {}/trimmed/reverse_reads.fasta.gz").format(self.concat_files['r'], self.outdir))

    batchfile.write("spades.py --threads {} --careful --memory {} -o {}/assembly"\
    .format(self.config["slurm_header"]["threads"], 8*int(self.config["slurm_header"]["threads"]), self.outdir))
    
    batchfile.write(" -1 {}/trimmed/forward_reads.fasta.gz".format(self.outdir))
    batchfile.write(" -2 {}/trimmed/reverse_reads.fasta.gz".format(self.outdir))
    batchfile.write(" -s {}".format(self.concat_files['i']))
    batchfile.write("rm {}/trimmed/forward_reads.fasta.gz {}/trimmed/reverse_reads.fasta.gz".format(self.outdir, self.outdir))

    batchfile.write("\n\n")
    batchfile.close()

  def create_resistancesection(self):
    """Creates a blast job for instances where many loci definition files make up an organism"""

    #Create run
    batchfile = open(self.batchfile, "a+")
    batchfile.write("mkdir {}/resistance\n\n".format(self.outdir))
    blast_format = "\"7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send length\""
    res_list = glob.glob("{}/*.fsa".format(self.config["folders"]["resistances"]))
    for entry in res_list:
      batchfile.write("# BLAST Resistance search in {} for {}\n".format(self.organism, os.path.basename(entry[:-4])))
      batchfile.write("blastn -db {}  -query {}/assembly/contigs.fasta -out {}/resistance/{}.txt -task megablast -num_threads {} -outfmt {}\n".format(\
      entry[:-4], self.outdir, self.outdir, os.path.basename(entry[:-4]), self.config["slurm_header"]["threads"], blast_format))
    batchfile.write("\n")
    batchfile.close()

  def create_mlstsection(self):
    """Creates a blast job for instances where many loci definition files make up an organism"""
    
    #Create run
    batchfile = open(self.batchfile, "a+")
    batchfile.write("mkdir {}/blast\n\n".format(self.outdir))
    blast_format = "\"7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send length\""
    tfa_list = glob.glob("{}/{}/*.tfa".format(self.config["folders"]["references"], self.organism))
    for entry in tfa_list:
      batchfile.write("# BLAST MLST alignment for {}, {}\n".format(self.organism, os.path.basename(entry[:-4])))
      batchfile.write("blastn -db {}  -query {}/assembly/contigs.fasta -out {}/blast/loci_query_{}.txt -task megablast -num_threads {} -outfmt {}\n"\
                      .format(entry[:-4], self.outdir, self.outdir, os.path.basename(entry[:-4]), self.config["slurm_header"]["threads"], blast_format))
    batchfile.write("\n")
    batchfile.close()

  def create_variantsection(self, trimmed=True):
    """ Creates a job for variant calling based on local alignment """
    ref = "{}/{}.fasta".format(self.config['folders']['genomes'],self.lims_fetcher.data['reference'])
    localdir = "{}/alignment".format(self.outdir)
    outbase = "{}/{}_{}".format(localdir, self.name, self.lims_fetcher.data['reference'])
    files = self.verify_fastq()

    if trimmed:
      reads_forward = self.concat_files['f'] 
      reads_reverse = self.concat_files['r']
    #State for old-MWGS comparison. Not normally reachable
    elif not trimmed:
      forward = list()
      reverse = list()
      for file in files[self.name]:
        fullfile = "{}/{}".format(self.indir, file)
        #Even indexes = Forward
        if not files[self.name].index(file)  % 2:
          forward.append(fullfile)
        elif files[self.name].index(file)  % 2:
          reverse.append(fullfile)
      reads_forward = "<( cat {} )".format(' '.join(forward)) 
      reads_reverse = "<( cat {} )".format(' '.join(reverse))
 
    #Create run
    batchfile = open(self.batchfile, "a+")
    batchfile.write("# Variant calling based on local alignment\n")
    batchfile.write("mkdir {}\n".format(localdir))

    batchfile.write("## Alignment & Deduplication\n")
    batchfile.write("bwa mem -M -t {} {} {} {} > {}.sam\n".format(self.config["slurm_header"]["threads"], ref ,reads_forward, reads_reverse, outbase))
    batchfile.write("samtools view --threads {} -b -o {}.bam -T {} {}.sam\n".format(self.config["slurm_header"]["threads"], outbase, ref, outbase))
    batchfile.write("samtools sort --threads {} -o {}.bam_sort {}.bam\n".format(self.config["slurm_header"]["threads"], outbase, outbase))
    batchfile.write("picard MarkDuplicates I={}.bam_sort O={}.bam_sort_rmdup M={}.stats.dup REMOVE_DUPLICATES=true\n".format(outbase, outbase, outbase))
    batchfile.write("samtools index {}.bam_sort_rmdup\n".format(outbase))
    batchfile.write("samtools idxstats {}.bam_sort_rmdup &> {}.stats.ref\n".format(outbase, outbase))
    #Removal of temp aligment files
    batchfile.write("rm {}.sam".format(outbase))
    batchfile.write("rm {}.bam".format(outbase))

    #Samtools duplicate calling, legacy
    #batchfile.write("samtools fixmate --threads {} -r -m {}.bam_sort {}.bam_sort_ms\n".format(self.config["slurm_header"]["threads"], outbase, outbase))
    #batchfile.write("samtools sort --threads {} -o {}.bam_sort {}.bam_sort_ms\n".format(self.config["slurm_header"]["threads"], outbase, outbase))
    #batchfile.write("samtools markdup -r -s --threads {} --reference {} --output-fmt bam {}.bam_sort {}.bam_sort_mkdup &> {}.stats.dup\n"\
    #                .format(self.config["slurm_header"]["threads"], ref, outbase, outbase, outbase))
    #batchfile.write("samtools rmdup --reference {} {}.bam_sort_mkdup {}.bam_sort_rmdup\n".format(ref, outbase, outbase))
    #batchfile.write("## Indexing\n")
    #batchfile.write("samtools index {}.bam_sort_mkdup\n".format(outbase))
    #batchfile.write("samtools idxstats {}.bam_sort_mkdup &> {}.stats.ref\n".format(outbase, outbase))

    batchfile.write("## Primary stats generation\n")
    #Insert stats, dedupped
    batchfile.write("samtools stats {}.bam_sort_rmdup |grep ^IS | cut -f 2- &> {}.stats.ins\n".format(outbase, outbase))
    #Coverage
    batchfile.write("samtools stats --coverage 1,10000,1 {}.bam_sort_rmdup |grep ^COV | cut -f 2- &> {}.stats.cov\n".format(outbase, outbase))
    #Mapped rate, no dedup,dedup in MWGS (trimming has no effect)!
    batchfile.write("samtools flagstat {}.bam_sort &> {}.stats.map\n".format(outbase, outbase))
    #Total reads, no dedup,dedup in MWGS (trimming has no effect)!
    batchfile.write("samtools view -c {}.bam_sort &> {}.stats.raw\n".format(outbase, outbase))

    batchfile.write("\n")
    batchfile.close()

  def create_trimsection(self):
    for root, dirs, files in os.walk(self.config["folders"]["adapters"]):
      if not "NexteraPE-PE.fa" in files:
        self.logger.error("Adapters folder at {} does not contain NexteraPE-PE.fa. Review paths.yml")
      else:
        break
    trimdir = "{}/trimmed".format(self.outdir)
    files = self.verify_fastq()
    batchfile = open(self.batchfile, "a+")
    batchfile.write("mkdir {}\n\n".format(trimdir))

    for index, k in enumerate(files):
      outfile = files[k][0].split('.')[0][:-2]
      if not outfile in self.trimmed_files:
        self.trimmed_files[outfile] = dict()
      self.trimmed_files[outfile]['fp'] = "{}/{}_trim_front_pair.fq".format(trimdir, outfile)
      self.trimmed_files[outfile]['fu'] = "{}/{}_trim_front_unpair.fq".format(trimdir, outfile)
      self.trimmed_files[outfile]['rp'] = "{}/{}_trim_rev_pair.fq".format(trimdir, outfile)
      self.trimmed_files[outfile]['ru'] = "{}/{}_trim_rev_unpair.fq".format(trimdir, outfile)

      batchfile.write("# Trimmomatic set {}\n".format(index+1))
      batchfile.write("trimmomatic PE -threads {} -phred33 {}/{} {}/{} {} {} {} {}\
      ILLUMINACLIP:{}NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\n\n"\
      .format(self.config["slurm_header"]["threads"], self.indir, files[k][0], self.indir, files[k][1],\
      self.trimmed_files[outfile]['fp'], self.trimmed_files[outfile]['fu'], self.trimmed_files[outfile]['rp'],\
      self.trimmed_files[outfile]['ru'], self.config["folders"]["adapters"]))

  def create_assemblystats_section(self):
    batchfile = open(self.batchfile, "a+")
    batchfile.write("# QUAST QC metrics\n")
    batchfile.write("mkdir {}/quast\n".format(self.outdir))
    batchfile.write("quast.py {}/assembly/contigs.fasta -o {}/quast\n\n".format(self.outdir, self.outdir))
    batchfile.close()

  def create_snpsection(self):
    snplist = self.filelist.copy()
    batchfile = open(self.batchfile, "a+")

    #VCFTools filters:
    vcffilter="--minQ 30 --thin 50 --minDP 3 --min-meanDP 20"
    #BCFTools filters:
    bcffilter = "GL[0]<-500 & GL[1]=0 & QR/RO>30 & QA/AO>30 & QUAL>5000 & ODDS>1100 & GQ>140 & DP>100 & MQM>59 & SAP<15 & PAIRED>0.9 & EPP>3"

    for item in snplist:
      name = item.split('/')[-2]
      if '_' in name:
        name = name.split('_')[0]
      self.lims_fetcher.load_lims_sample_info(name)
      batchfile.write('# Basecalling for sample {}\n'.format(name))
      ref = "{}/{}.fasta".format(self.config['folders']['genomes'],self.lims_fetcher.data['reference'])
      outbase = "{}/{}_{}".format(item, name, self.lims_fetcher.data['reference'])
      batchfile.write("samtools view -h -q 1 -F 4 -F 256 {}.bam_sort_rmdup | grep -v XA:Z | grep -v SA:Z| samtools view -b - > {}/{}.unique\n".format(outbase, self.outdir, name))
      batchfile.write('freebayes -= --pvar 0.7 -j -J --standard-filters -C 6 --min-coverage 30 --ploidy 1 -f {} -b {}/{}.unique -v {}/{}.vcf\n'.format(ref, self.outdir, name , self.outdir, name))
      batchfile.write('bcftools view {}/{}.vcf -o {}/{}.bcf.gz -O b --exclude-uncalled --types snps\n'.format(self.outdir, name, self.outdir, name))
      batchfile.write('bcftools index {}/{}.bcf.gz\n'.format(self.outdir, name))
      batchfile.write('\n')

      batchfile.write('vcftools --bcf {}/{}.bcf.gz {} --remove-filtered-all --recode-INFO-all --recode-bcf --out {}/{}\n'.format(self.outdir, name, vcffilter, self.outdir, name))
      batchfile.write('bcftools view {}/{}.recode.bcf -i "{}" -o {}/{}.recode.bcf.gz -O b --exclude-uncalled --types snps\n'.format(self.outdir, name, bcffilter, self.outdir, name))
      batchfile.write('bcftools index {}/{}.recode.bcf.gz\n\n'.format(self.outdir, name))

    batchfile.write('# SNP pair-wise distance\n')
    batchfile.write('touch {}/stats.out\n'.format(self.outdir))
    while len(snplist) > 1:
      top = snplist.pop(0)
      nameOne = top.split('/')[-2]
      if '_' in nameOne:
        nameOne = nameOne.split('_')[0]
      for entry in snplist:
        nameTwo = entry.split('/')[-2]
        if '_' in nameTwo:
          nameTwo = nameTwo.split('_')[0]

        pair = "{}_{}".format(nameOne, nameTwo)
        batchfile.write('bcftools isec {}/{}.recode.bcf.gz {}/{}.recode.bcf.gz -n=1 -c all -p {}/tmp -O b\n'.format(self.outdir, nameOne, self.outdir, nameTwo, self.outdir))
        batchfile.write('bcftools merge -O b -o {}/{}.bcf.gz --force-samples {}/tmp/0000.bcf {}/tmp/0001.bcf\n'.format(self.outdir, pair, self.outdir, self.outdir))
        batchfile.write('bcftools index {}/{}.bcf.gz\n'.format(self.outdir, pair))

        batchfile.write("echo {} $( bcftools stats {}/{}.bcf.gz |grep SNPs: | cut -d $'\\t' -f4 ) >> {}/stats.out\n".format(pair, self.outdir, pair, self.outdir))
        batchfile.write('\n')
    batchfile.close()


  def create_project(self, name):
    """Creates project in database"""
    try:
      self.lims_fetcher.load_lims_project_info(name)
    except Exception as e:
      self.logger.error("Unable to load LIMS info for project {}".format(name))
    proj_col=dict()
    proj_col['CG_ID_project'] = name
    proj_col['Customer_ID_project'] = self.lims_fetcher.data['Customer_ID_project']
    proj_col['date_ordered'] = self.lims_fetcher.data['date_received']
    self.db_pusher.add_rec(proj_col, 'Projects')

  def create_sample(self, name):
    """Creates sample in database"""
    try:
      self.lims_fetcher.load_lims_sample_info(name)
      sample_col = self.db_pusher.get_columns('Samples') 
      sample_col['CG_ID_sample'] = self.lims_fetcher.data['CG_ID_sample']
      sample_col['CG_ID_project'] = self.lims_fetcher.data['CG_ID_project']
      sample_col['Customer_ID_sample'] = self.lims_fetcher.data['Customer_ID_sample']
      sample_col['reference_genome'] = self.lims_fetcher.data['reference']
      sample_col["date_analysis"] = self.dt
      sample_col['organism']=self.lims_fetcher.data['organism']
      #self.db_pusher.purge_rec(sample_col['CG_ID_sample'], 'sample')
      self.db_pusher.add_rec(sample_col, 'Samples')
    except Exception as e:
      self.logger.error("Unable to add sample {} to database".format(self.name))

  def project_job(self, single_sample=False, qc_only=False, trimmed=True):
    if 'dry' in self.config and self.config['dry']==True:
      dry=True
    else:
      dry=False
    jobarray = list()
    if not os.path.exists(self.finishdir):
      os.makedirs(self.finishdir)
    try:
       if single_sample:
         self.create_project(os.path.normpath(self.indir).split('/')[-2])
       else:
        self.create_project(self.name)
    except Exception as e:
      self.logger.error("LIMS interaction failed. Unable to read/write project {}".format(self.name))
      #Start every sample job
    if single_sample:
      try:
        self.sample_job(qc_only=qc_only, trimmed=trimmed)
        headerargs = self.get_headerargs()
        outfile = self.get_sbatch()
        bash_cmd="sbatch {} {}".format(headerargs, outfile)
        if not dry and outfile != "":
          samproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
          output, error = samproc.communicate()
          jobno = re.search('(\d+)', str(output)).group(0)
          jobarray.append(jobno)
        else:
          self.logger.info("Suppressed command: {}".format(bash_cmd))
      except Exception as e:
        self.logger.error("Unable to analyze single sample {}".format(self.name))
    else:
      for (dirpath, dirnames, filenames) in os.walk(self.indir):
        for dir in dirnames:
          try:
            sample_in = "{}/{}".format(dirpath, dir)
            sample_out = "{}/{}".format(self.finishdir, dir)
            sample_instance = Job_Creator(sample_in, self.config, self.logger, sample_out, self.now) 
            sample_instance.sample_job(qc_only=qc_only, trimmed=trimmed)
            headerargs = sample_instance.get_headerargs()
            outfile = ""
            if os.path.isfile(sample_instance.get_sbatch()):
              outfile = sample_instance.get_sbatch()
            bash_cmd="sbatch {} {}".format(headerargs, outfile)
            if not dry and outfile != "":
              projproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
              output, error = projproc.communicate()
              jobno = re.search('(\d+)', str(output)).group(0)
              jobarray.append(jobno)
            else:
              self.logger.info("Suppressed command: {}".format(bash_cmd))
          except Exception as e:
            pass
    if not dry:
      self.finish_job(jobarray, single_sample)

  def finish_job(self, joblist, single_sample=False):
    """ Uploads data and sends an email once all analysis jobs are complete. """

    startfile = "{}/run_started.out".format(self.finishdir)
    mailfile = "{}/mailjob.sh".format(self.finishdir)
    mb = open(mailfile, "w+")
    sb = open(startfile, "w+")
    sb.write("#!/usr/bin/env bash\n\n")
    sb.close()
    mb.write("#!/usr/bin/env bash\n\n")
    mb.write("#Uploading of results to database and production of report\n")
    if 'MICROSALT_CONFIG' in os.environ:
      mb.write("export MICROSALT_CONFIG={}\n".format(os.environ['MICROSALT_CONFIG']))
    mb.write("source activate $CONDA_DEFAULT_ENV\n")
    if not single_sample:
      mb.write("microSALT utils finish project {} --input {} --rerun --email {}\n".\
               format(self.name, self.finishdir, self.config['regex']['mail_recipient']))
    else:
      mb.write("microSALT utils finish sample {} --input {} --rerun --email {}\n".\
               format(self.name, self.finishdir, self.config['regex']['mail_recipient']))
    mb.write("touch {}/run_complete.out".format(self.finishdir))
    mb.close()

    massagedJobs = list()
    final = ':'.join(joblist)
    #Create subtracker if more than 50 samples
    maxlen = 50
    if len(joblist) > maxlen:
      i = 1
      while i <= len(joblist):
        if i+maxlen < len(joblist):
          massagedJobs.append(':'.join(joblist[i-1:i+maxlen-1]))
        else:
          massagedJobs.append(':'.join(joblist[i-1:-1]))
        i += maxlen
      for entry in massagedJobs:
        if massagedJobs.index(entry) < len(massagedJobs)-1:
          head = "-A {} -p core -n 1 -t 00:00:10 -J {}_{}_SUBTRACKER --qos {} --dependency=afterany:{}"\
                 .format(self.config["slurm_header"]["project"],self.config["slurm_header"]["job_prefix"],\
                         self.name,self.config["slurm_header"]["qos"],entry)
          bash_cmd="sbatch {} {}".format(head, startfile)
          mailproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
          output, error = mailproc.communicate()
          jobno = re.search('(\d+)', str(output)).group(0)
          massagedJobs[massagedJobs.index(entry)+1] += ":{}".format(jobno)
        else:
          final = entry
          break 

    head = "-A {} -p core -n 1 -t 06:00:00 -J {}_{}_MAILJOB --qos {} --open-mode append --dependency=afterany:{} --output {}"\
            .format(self.config["slurm_header"]["project"],self.config["slurm_header"]["job_prefix"],\
                    self.name,self.config["slurm_header"]["qos"],\
           final, self.config['folders']['log_file'],  self.config['regex']['mail_recipient'])
    bash_cmd="sbatch {} {}".format(head, mailfile)
    mailproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
    output, error = mailproc.communicate()

  def sample_job(self, qc_only=False, trimmed=True):
    """ Writes necessary sbatch job for each individual sample """
    try:
      self.trimmed_files = dict()
      if not os.path.exists(self.finishdir):
        os.makedirs(self.finishdir)
      try:
        self.organism = self.lims_fetcher.get_organism_refname(self.name, external=False)
        # This is one job 
        self.batchfile = "{}/runfile.sbatch".format(self.finishdir)
        batchfile = open(self.batchfile, "w+")
        batchfile.write("#!/usr/bin/env bash\n\n")
        batchfile.write("mkdir -p {}\n".format(self.outdir))
        batchfile.close()

        self.create_trimsection()
        self.interlace_files()
        #self.logger.info("Sample trimming is currently disabled for QC results")
        self.create_variantsection(trimmed=trimmed)
        if not qc_only:
          self.create_assemblysection()
          self.create_assemblystats_section()
          self.create_mlstsection()
          self.create_resistancesection()
        batchfile = open(self.batchfile, "a+")
        batchfile.write("cp -r {}/* {}".format(self.outdir, self.finishdir))
        batchfile.close()

        self.logger.info("Created runfile for sample {} in folder {}".format(self.name, self.outdir))
      except Exception as e:
        raise 
      try: 
        self.create_sample(self.name)
      except Exception as e:
        self.logger.error("Unable to access LIMS info for sample {}".format(self.name))
    except Exception as e:
      self.logger.error("Unable to create job for sample {}\nSource: {}".format(self.name, str(e)))
      shutil.rmtree(self.finishdir, ignore_errors=True)
      raise

  def snp_job(self):
    """ Writes a SNP calling job for a set of samples """
    if not os.path.exists(self.finishdir):
      os.makedirs(self.finishdir)

    self.batchfile = "{}/runfile.sbatch".format(self.finishdir)
    batchfile = open(self.batchfile, "w+")
    batchfile.write("#!/usr/bin/env bash\n\n")
    batchfile.write("mkdir -p {}\n".format(self.outdir))
    batchfile.close()

    self.create_snpsection()
    batchfile = open(self.batchfile, "a+")
    batchfile.write("cp -r {}/* {}".format(self.outdir, self.finishdir))
    batchfile.close()

    headerline = "-A {} -p {} -n 1 -t 24:00:00 -J {}_{} --qos {} --output {}/slurm_{}.log".format(self.config["slurm_header"]["project"],\
                 self.config["slurm_header"]["type"],\
                 self.config["slurm_header"]["job_prefix"], self.name,self.config["slurm_header"]["qos"], self.finishdir, self.name)
    outfile = self.get_sbatch()
    bash_cmd="sbatch {} {}".format(headerline, outfile)
    samproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
    output, error = samproc.communicate()
