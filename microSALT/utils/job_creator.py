"""Creates sbatch jobs for MLST instances
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import glob
import json
import os
import re
import shutil
import subprocess
import time

from datetime import datetime
from microSALT.store.lims_fetcher import LIMS_Fetcher
from microSALT.store.db_manipulator import DB_Manipulator

class Job_Creator():

  def __init__(self, input, config, log, finishdir="", timestamp="", trim=True, qc_only=False,careful=False,pool=list()):
    self.config = config
    self.logger = log
    self.batchfile = ""
    self.filelist = list()
    self.indir = ""
    self.trimmed=trim
    self.qc_only = qc_only
    self.careful = careful
    self.pool = pool

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

    self.finishdir = finishdir
    if self.finishdir == "":
      self.finishdir="{}/{}_{}".format(config["folders"]["results"], self.name, self.now)
    self.db_pusher=DB_Manipulator(config, log)
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
    verified_files = list()
    files = os.listdir(self.indir)
    if files == []:
      raise Exception("Directory {} lacks fastq files.".format(self.indir))
    for file in files:
      file_match = re.match( self.config['regex']['file_pattern'], file)
      if file_match:
        #Check that symlinks resolve
        path = '{}/{}'.format(self.indir, file)
        if os.path.islink(path):
          if not os.path.exists(os.readlink(path)):
            raise Exception("Some fastq files are unresolved symlinks in directory {}.".format(self.indir))

        #Make sure both mates exist
        if file_match[1] == '1' or file_match[1] == '2':
          pairno = 2-1%int(file_match[1]) #1->2, 2->1
          #Construct mate name
          pairname = "{}{}{}".format(file_match.string[:file_match.end(1)-1] , pairno, \
                      file_match.string[file_match.end(1):file_match.end()])
          if pairname in files:
            files.pop( files.index(pairname) )
            verified_files.append(file_match[0])
            verified_files.append(pairname)
        else:
          raise Exception("Some fastq files have no mate in directory {}.".format(self.indir))
    if verified_files == []:
      raise Exception("No files in directory {} match file_pattern '{}'.".format(self.indir, self.config['regex']['file_pattern']))
    #Warn about file sizes
    for vfile in verified_files:
      bsize = os.stat("{}/{}".format(self.indir,vfile)).st_size
      bsize = bsize >> 20
      if bsize > 1000:
        self.logger.warning("Input fastq {} exceeds 1000MB".format(vfile))
    return sorted(verified_files)
 
  def create_assemblysection(self):
    batchfile = open(self.batchfile, "a+")
    #memory is actually 128 per node regardless of cores.
    batchfile.write("# Spades assembly\n")
    if self.trimmed:
      trimline = '-s {}'.format(self.concat_files['i'])
    else:
      trimline = ''
    if self.careful:
      careline = '--careful'
    else:
      careline = ''
    
    batchfile.write("spades.py --threads {} {} --memory {} -o {}/assembly -1 {} -2 {} {}\n"\
    .format(self.config["slurm_header"]["threads"], careline, 8*int(self.config["slurm_header"]["threads"]), self.finishdir, self.concat_files['f'], self.concat_files['r'], trimline))
    batchfile.write("##Input cleanup\n")
    batchfile.write("rm -r {}/trimmed\n".format(self.finishdir))
    batchfile.write("\n\n")
    batchfile.close()

  def blast_subset(self, name, search_string):
    #Create run
    file_list = glob.glob(search_string)
    batchfile = open(self.batchfile, "a+")
    batchfile.write("mkdir {}/blast_search/{}\n".format(self.finishdir, name))
    header_format = "\"7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send length\""
    id_format = "\"7 sallseqid sstrand qaccver saccver pident evalue bitscore qstart qend sstart send length\""

    if name == 'cgmlst':
      blast_format = id_format
    else:
      blast_format = header_format

    if len(file_list) > 1:
      for ref in file_list:
        ref_nosuf = re.search('(\w+)\.\w+', os.path.basename(ref)).group(1)
        batchfile.write("# BLAST {} search for {}, {}\n".format(name, self.organism, ref_nosuf))
        if name == 'mlst':
          batchfile.write("blastn -db {}/{}  -query {}/assembly/contigs.fasta -out {}/blast_search/{}/loci_query_{}.txt -task megablast -num_threads {} -outfmt {}\n".format(\
          os.path.dirname(ref), ref_nosuf, self.finishdir, self.finishdir, name, ref_nosuf, self.config["slurm_header"]["threads"], blast_format))
        else:
          batchfile.write("blastn -db {}/{}  -query {}/assembly/contigs.fasta -out {}/blast_search/{}/{}.txt -task megablast -num_threads {} -outfmt {}\n".format(\
          os.path.dirname(ref), ref_nosuf, self.finishdir, self.finishdir, name, ref_nosuf, self.config["slurm_header"]["threads"], blast_format))
    else:
      ref_nosuf = re.search('(\w+)\.\w+', os.path.basename(file_list[0])).group(1)
      batchfile.write("## BLAST {} search in {}\n".format(name, self.organism.replace('_', ' ').capitalize() ))
      batchfile.write("blastn -db {}/{}  -query {}/assembly/contigs.fasta -out {}/blast_search/{}/{}.txt -task megablast -num_threads {} -outfmt {}\n".format(\
                    os.path.dirname(search_string), ref_nosuf, self.finishdir, self.finishdir, name, ref_nosuf, self.config["slurm_header"]["threads"], blast_format))
    batchfile.write("\n")
    batchfile.close()

  def create_variantsection(self):
    """ Creates a job for variant calling based on local alignment """
    ref = "{}/{}.fasta".format(self.config['folders']['genomes'],self.lims_fetcher.data['reference'])
    localdir = "{}/alignment".format(self.finishdir)
    outbase = "{}/{}_{}".format(localdir, self.name, self.lims_fetcher.data['reference'])

    #Create run
    batchfile = open(self.batchfile, "a+")
    batchfile.write("# Variant calling based on local alignment\n")
    batchfile.write("mkdir {}\n".format(localdir))

    batchfile.write("## Alignment & Deduplication\n")
    batchfile.write("bwa mem -M -t {} {} {} {} > {}.sam\n".format(self.config["slurm_header"]["threads"], ref ,self.concat_files['f'], self.concat_files['r'], outbase))
    batchfile.write("samtools view --threads {} -b -o {}.bam -T {} {}.sam\n".format(self.config["slurm_header"]["threads"], outbase, ref, outbase))
    batchfile.write("samtools sort --threads {} -o {}.bam_sort {}.bam\n".format(self.config["slurm_header"]["threads"], outbase, outbase))
    batchfile.write("picard MarkDuplicates I={}.bam_sort O={}.bam_sort_rmdup M={}.stats.dup REMOVE_DUPLICATES=true\n".format(outbase, outbase, outbase))
    batchfile.write("samtools index {}.bam_sort_rmdup\n".format(outbase))
    batchfile.write("samtools idxstats {}.bam_sort_rmdup &> {}.stats.ref\n".format(outbase, outbase))
    #Removal of temp aligment files
    batchfile.write("rm {}.bam {}.sam\n".format(outbase, outbase))

    batchfile.write("## Primary stats generation\n")
    #Insert stats, dedupped
    batchfile.write("picard CollectInsertSizeMetrics I={}.bam_sort_rmdup O={}.stats.ins H={}.hist.ins\n".format(outbase, outbase, outbase))
    #Coverage
    batchfile.write("samtools stats --coverage 1,10000,1 {}.bam_sort_rmdup |grep ^COV | cut -f 2- &> {}.stats.cov\n".format(outbase, outbase))
    #Mapped rate, no dedup,dedup in MWGS (trimming has no effect)!
    batchfile.write("samtools flagstat {}.bam_sort &> {}.stats.map\n".format(outbase, outbase))
    #Total reads, no dedup,dedup in MWGS (trimming has no effect)!
    batchfile.write("samtools view -c {}.bam_sort &> {}.stats.raw\n".format(outbase, outbase))

    batchfile.write("\n\n")
    batchfile.close()

  def create_preprocsection(self):
    """Concatinates data, possibly trims it, then makes the unstranded reads usable"""
    forward = list()
    reverse = list()

    for root, dirs, files in os.walk(self.config["folders"]["adapters"]):
      if not "NexteraPE-PE.fa" in files:
        self.logger.error("Adapters folder at {} does not contain NexteraPE-PE.fa. Review paths.yml")
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
      #Even indexes = Forward
      if not files.index(file)  % 2:
        forward.append(fullfile)
      elif files.index(file)  % 2:
        reverse.append(fullfile)
    outfile = files[0].split('_')[0]

    self.concat_files['f'] = "{}/trimmed/forward_reads.fastq.gz".format(self.finishdir)
    self.concat_files['r'] = "{}/trimmed/reverse_reads.fastq.gz".format(self.finishdir)
    batchfile.write("cat {} > {}\n".format(' '.join(forward), self.concat_files['f']))
    batchfile.write("cat {} > {}\n".format(' '.join(reverse), self.concat_files['r']))

    if self.trimmed:
      fp = "{}/{}_trim_front_pair.fastq.gz".format(trimdir, outfile)
      fu = "{}/{}_trim_front_unpair.fastq.gz".format(trimdir, outfile)
      rp = "{}/{}_trim_rev_pair.fastq.gz".format(trimdir, outfile)
      ru = "{}/{}_trim_rev_unpair.fastq.gz".format(trimdir, outfile)
      batchfile.write("##Trimming section\n")
      batchfile.write("trimmomatic PE -threads {} -phred33 {} {} {} {} {} {}\
      ILLUMINACLIP:{}/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\n"\
      .format(self.config["slurm_header"]["threads"], self.concat_files['f'], self.concat_files['r'], fp, fu, rp, ru, self.config["folders"]["adapters"]))

      batchfile.write("## Interlaced trimmed files\n")
      self.concat_files['f'] = fp
      self.concat_files['r'] = rp 
      self.concat_files['i'] = "{}/{}_trim_unpair.fastq.gz".format(trimdir, outfile)

      batchfile.write("cat {} >> {}\n".format(' '.join([fu, ru]), self.concat_files['i']))
    batchfile.write("\n")
    batchfile.close()

  def create_assemblystats_section(self):
    batchfile = open(self.batchfile, "a+")
    batchfile.write("# QUAST QC metrics\n")
    batchfile.write("mkdir {}/assembly/quast\n".format(self.finishdir))
    batchfile.write("quast.py {}/assembly/contigs.fasta -o {}/assembly/quast\n\n".format(self.finishdir, self.finishdir))
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
      batchfile.write("samtools view -h -q 1 -F 4 -F 256 {}.bam_sort_rmdup | grep -v XA:Z | grep -v SA:Z| samtools view -b - > {}/{}.unique\n".format(outbase, self.finishdir, name))
      batchfile.write('freebayes -= --pvar 0.7 -j -J --standard-filters -C 6 --min-coverage 30 --ploidy 1 -f {} -b {}/{}.unique -v {}/{}.vcf\n'.format(ref, self.finishdir, name , self.finishdir, name))
      batchfile.write('bcftools view {}/{}.vcf -o {}/{}.bcf.gz -O b --exclude-uncalled --types snps\n'.format(self.finishdir, name, self.finishdir, name))
      batchfile.write('bcftools index {}/{}.bcf.gz\n'.format(self.finishdir, name))
      batchfile.write('\n')

      batchfile.write('vcftools --bcf {}/{}.bcf.gz {} --remove-filtered-all --recode-INFO-all --recode-bcf --out {}/{}\n'.format(self.finishdir, name, vcffilter, self.finishdir, name))
      batchfile.write('bcftools view {}/{}.recode.bcf -i "{}" -o {}/{}.recode.bcf.gz -O b --exclude-uncalled --types snps\n'.format(self.finishdir, name, bcffilter, self.finishdir, name))
      batchfile.write('bcftools index {}/{}.recode.bcf.gz\n\n'.format(self.finishdir, name))

    batchfile.write('# SNP pair-wise distance\n')
    batchfile.write('touch {}/stats.out\n'.format(self.finishdir))
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
        batchfile.write('bcftools isec {}/{}.recode.bcf.gz {}/{}.recode.bcf.gz -n=1 -c all -p {}/tmp -O b\n'.format(self.finishdir, nameOne, self.finishdir, nameTwo, self.finishdir))
        batchfile.write('bcftools merge -O b -o {}/{}.bcf.gz --force-samples {}/tmp/0000.bcf {}/tmp/0001.bcf\n'.format(self.finishdir, pair, self.finishdir, self.finishdir))
        batchfile.write('bcftools index {}/{}.bcf.gz\n'.format(self.finishdir, pair))

        batchfile.write("echo {} $( bcftools stats {}/{}.bcf.gz |grep SNPs: | cut -d $'\\t' -f4 ) >> {}/stats.out\n".format(pair, self.finishdir, pair, self.finishdir))
        batchfile.write('\n')
    batchfile.close()

  def create_collection(self):
    """Creates collection entry in database"""
    if self.db_pusher.exists('Collections', {'ID_collection':self.name}):
      self.db_pusher.purge_rec(name=self.name, type='Collections')
      for sample in self.pool:
        self.db_pusher.add_rec({'ID_collection':self.name, 'CG_ID_sample':sample}, 'Collections')

    addedprojs = list()
    for sample in self.pool:
      proj = re.search('(\w+)A(?:\w+)', sample).group(1)
      if proj not in addedprojs:
        self.create_project(proj)
        addedprojs.append(proj)

  def create_project(self, name):
    """Creates project in database"""
    try:
      self.lims_fetcher.load_lims_project_info(name)
    except Exception as e:
      self.logger.error("Unable to load LIMS info for project {}".format(name))
    proj_col=dict()
    proj_col['CG_ID_project'] = name
    proj_col['Customer_ID_project'] = self.lims_fetcher.data['Customer_ID_project']
    proj_col['Customer_ID'] = self.lims_fetcher.data['Customer_ID']
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
      sample_col["application_tag"] = self.lims_fetcher.data['application_tag']
      sample_col["priority"] = self.lims_fetcher.data['priority']
      sample_col["date_arrival"] = self.lims_fetcher.data['date_arrival']
      sample_col["date_sequencing"] = self.lims_fetcher.data['date_sequencing']
      sample_col["date_libprep"] = self.lims_fetcher.data['date_libprep']
      sample_col["method_libprep"] = self.lims_fetcher.data['method_libprep']
      sample_col["method_sequencing"] = self.lims_fetcher.data['method_sequencing']
      #self.db_pusher.purge_rec(sample_col['CG_ID_sample'], 'sample')
      self.db_pusher.add_rec(sample_col, 'Samples')
    except Exception as e:
      self.logger.error("Unable to add sample {} to database".format(self.name))

  def project_job(self, single_sample=False):
    if 'dry' in self.config and self.config['dry']==True:
      dry=True
    else:
      dry=False
    jobarray = list()
    if not os.path.exists(self.finishdir):
      os.makedirs(self.finishdir)
    #Loads project level info.
    try:
       if single_sample:
         self.create_project(os.path.normpath(self.indir).split('/')[-2])
       elif self.pool:
         self.create_collection()
       else:
        self.create_project(self.name)
    except Exception as e:
      self.logger.error("LIMS interaction failed. Unable to read/write project {}".format(self.name))
    #Writes the job creation sbatch
    if single_sample:
      try:
        self.sample_job()
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
    elif self.pool:
      for (dirpath, dirnames, filenames) in os.walk(self.indir):
        for dir in dirnames:
          try:
            sample_in = "{}/{}".format(dirpath, dir)
            sample_out = "{}/{}".format(self.finishdir, dir)
            sample_instance = Job_Creator(sample_in, self.config, self.logger, sample_out, self.now, trim=self.trimmed, careful=self.careful)
            sample_instance.sample_job()
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
    else:
      for (dirpath, dirnames, filenames) in os.walk(self.indir):
        for dir in dirnames:
          try:
            sample_in = "{}/{}".format(dirpath, dir)
            sample_out = "{}/{}".format(self.finishdir, dir)
            sample_instance = Job_Creator(sample_in, self.config, self.logger, sample_out, self.now, trim=self.trimmed, careful=self.careful) 
            sample_instance.sample_job()
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
    report = 'default'
    scope = 'project'
    extraflag = '--rerun'
    if single_sample:
      scope = 'sample' 
    elif self.pool:
      scope = 'collection'
      report = 'motif_overview'
      extraflag = ''
    if self.qc_only:
      report = 'qc'
    custom_conf = ''
    if 'config_path' in self.config:
      custom_conf = '--config {}'.format(self.config['config_path'])

    process = subprocess.Popen("id -un".split(), stdout=subprocess.PIPE)
    user, error = process.communicate()
    user = str(user).replace('.',' ').title()
    startfile = "{}/run_started.out".format(self.finishdir)
    configfile = "{}/config.log".format(self.finishdir) 
    mailfile = "{}/mailjob.sh".format(self.finishdir)
    sb = open(startfile, "w+")
    cb = open(configfile, "w+")
    mb = open(mailfile, "w+")
    sb.write("#!/usr/bin/env bash\n")
    sb.close()
    configout = self.config.copy()
    if 'genologics' in configout:
      del configout['genologics']
    cb.write("ANALYSIS STARTED BY: {}\n".format(user))
    cb.write(json.dumps(configout, indent=2,separators=(',',':')))
    cb.close()
    mb.write("#!/usr/bin/env bash\n\n")
    mb.write("#Uploading of results to database and production of report\n")
    if 'MICROSALT_CONFIG' in os.environ:
      mb.write("export MICROSALT_CONFIG={}\n".format(os.environ['MICROSALT_CONFIG']))
    mb.write("source activate $CONDA_DEFAULT_ENV\n")

    mb.write("microSALT utils finish {} {} --input {} {} --email {} --report {} {}\n".\
               format(scope, self.name, self.finishdir, extraflag, self.config['regex']['mail_recipient'], report, custom_conf))
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

    head = "-A {} -p core -n 1 -t 18:00:00 -J {}_{}_MAILJOB --qos {} --open-mode append --dependency=afterany:{} --output {}"\
            .format(self.config["slurm_header"]["project"],self.config["slurm_header"]["job_prefix"],\
                    self.name,self.config["slurm_header"]["qos"],\
           final, self.config['folders']['log_file'],  self.config['regex']['mail_recipient'])
    bash_cmd="sbatch {} {}".format(head, mailfile)
    mailproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
    output, error = mailproc.communicate()

  def sample_job(self):
    """ Writes necessary sbatch job for each individual sample """
    try:
      if not os.path.exists(self.finishdir):
        os.makedirs(self.finishdir)
      try:
        self.organism = self.lims_fetcher.get_organism_refname(self.name)
        if not self.organism:
          self.organism = self.lims_fetcher.data['organism']
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

        self.logger.info("Created runfile for sample {} in folder {}".format(self.name, self.finishdir))
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

  def create_blast_search(self):
    self.batchfile = "{}/runfile.sbatch".format(self.finishdir)
    batchfile = open(self.batchfile, "a+")
    batchfile.write("mkdir -p {}/blast_search\n".format(self.finishdir))
    batchfile.close()
    self.blast_subset('mlst',"{}/{}/*.tfa".format(self.config["folders"]["st_loci"], self.organism))
    self.blast_subset('resistance',"{}/*.fsa".format(self.config["folders"]["resistances"]))
    if os.path.exists("{}/{}".format(self.config['folders']['cgmlst'], self.organism)):
      self.blast_subset('cgmlst', "{}/{}/*.fasta".format(self.config['folders']['cgmlst'], self.organism))
    self.blast_subset('fimh',"{}/*.fsa".format(self.config["folders"]["fimhs"]))
    self.blast_subset('plasmid',"{}/*.fsa".format(self.config["folders"]["plasmids"]))
    self.blast_subset('virulence',"{}/*.fsa".format(self.config["folders"]["virulences"]))
    if self.organism == "escherichia_coli":
      self.blast_subset('expec', self.config["folders"]["expec"])

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

    headerline = "-A {} -p {} -n 1 -t 24:00:00 -J {}_{} --qos {} --output {}/slurm_{}.log".format(self.config["slurm_header"]["project"],\
                 self.config["slurm_header"]["type"],\
                 self.config["slurm_header"]["job_prefix"], self.name,self.config["slurm_header"]["qos"], self.finishdir, self.name)
    outfile = self.get_sbatch()
    bash_cmd="sbatch {} {}".format(headerline, outfile)
    samproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
    output, error = samproc.communicate()
