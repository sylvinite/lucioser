"""Scrapes output files for data and adds them to the database
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import glob
import os
import re
import string
import sys
import time

from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.store.lims_fetcher import LIMS_Fetcher
from microSALT.utils.reporter import Reporter
from microSALT.utils.job_creator import Job_Creator

# TODO: Rewrite so samples use seperate objects
class Scraper():

  def __init__(self, infolder, config, log):
    self.config = config
    self.logger = log
    self.db_pusher=DB_Manipulator(config, log)
    self.lims_fetcher=LIMS_Fetcher(config, log)
    self.job_fallback=Job_Creator("", config, log)
    self.infolder = os.path.abspath(infolder)
    self.sampledir = ""
   
    last_folder = self.infolder.split('/')[-1]
    self.name = last_folder.split('_')[0]
    #Cruddy internal check
    if not 'ACC' in self.name and not 'MIC' in self.name:
      self.lims_fetcher.load_lims_sample_info(self.name)
      self.name = self.lims_fetcher.data['CG_ID_sample']
    #TODO: Replace date from dir with entry from analysis files/database
    if not '_' in last_folder:
      last_folder = self.infolder.split('/')[-2]
    self.date = "{} {}".format(re.sub(r'\.','-', last_folder.split('_')[1]), re.sub(r'\.',':', last_folder.split('_')[2]))
    self.lims_sample_info = {}
    self.gene2resistance = self.load_resistances()

  def scrape_project(self):
    """Scrapes a project folder for information"""
    if self.config['rerun']:
      self.db_pusher.purge_rec(self.name, 'Projects')
    if not self.db_pusher.exists('Projects', {'CG_ID_project':self.name}):
      self.logger.warning("Replacing project {}".format(self.name))
      self.job_fallback.create_project(self.name)

    #Scrape order matters a lot!
    for dir in os.listdir(self.infolder):
     if os.path.isdir("{}/{}".format(self.infolder, dir)): 
       self.sampledir = "{}/{}".format(self.infolder, dir)
       self.name = dir
       if not self.db_pusher.exists('Samples', {'CG_ID_sample':self.name}):
         self.logger.warning("Replacing sample {}".format(self.name))
         self.job_fallback.create_sample(self.name)
       self.scrape_blast(type='seq_type')
       self.scrape_blast(type='resistance')
       self.scrape_blast(type='core_seq_type')
       if self.lims_fetcher.get_organism_refname(self.name) == "escherichia_coli":
         self.scrape_blast(type='expec')
       self.scrape_alignment()
       self.scrape_quast()

  def scrape_sample(self):
    """Scrapes a sample folder for information"""
    if self.config['rerun']:
      self.db_pusher.purge_rec(self.name, 'Samples')

    self.lims_fetcher.load_lims_sample_info(self.name)
    if not self.db_pusher.exists('Projects', {'CG_ID_project':self.lims_fetcher.data['CG_ID_project']}):
      self.logger.warning("Replacing project {}".format(self.lims_fetcher.data['CG_ID_project']))
      self.job_fallback.create_project(self.lims_fetcher.data['CG_ID_project'])

    if not self.db_pusher.exists('Samples', {'CG_ID_sample':self.name}):
      self.logger.warning("Replacing sample {}".format(self.name))
      self.job_fallback.create_sample(self.name)

    #Scrape order matters a lot!
    self.sampledir = self.infolder
    self.scrape_blast(type='seq_type')
    self.scrape_blast(type='resistance')
    self.scrape_blast(type='core_seq_type')
    if self.lims_fetcher.get_organism_refname(self.name) == "escherichia_coli":
      self.scrape_blast(type='expec')
    self.scrape_alignment()
    self.scrape_quast()

  def scrape_quast(self):
    """Scrapes a quast report for assembly information"""
    quast = dict()
    report = "{}/assembly/quast/report.tsv".format(self.sampledir)
    try:
      with open(report, 'r') as infile:
        for line in infile:
          lsplit = line.rstrip().split('\t')
          if lsplit[0] == '# contigs':
            quast['contigs'] = int(lsplit[1])
          elif lsplit[0] == 'Total length':
            quast['genome_length'] = int(lsplit[1])
          elif lsplit[0] == 'GC (%)':
            quast['gc_percentage'] = float(lsplit[1])
          elif lsplit[0] == 'N50':
            quast['n50'] = int(lsplit[1])

      self.db_pusher.upd_rec({'CG_ID_sample' : self.name}, 'Samples', quast)
      # Too spammy
      #self.logger.info("Project {} recieved quast stats: {}"\
      #               .format(self.name, quast))
    except Exception as e:
      self.logger.warning("Cannot generate quast statistics for {}".format(self.name))

  def get_locilengths(self, filename):
    """ Generate a dict of length for any given loci """
    #Create dict with full name as key, associated nucleotides as value. 
    alleles=dict()
    finalalleles=dict()
    lastallele=""
    f = open(filename,"r")
    for row in f:
      if ">" in row:
        lastallele = row.strip()
        alleles[lastallele] = ""
      else:
        alleles[lastallele] = alleles[lastallele] + row.strip()
    f.close()

    for k, v in alleles.items():
      finalalleles[k] = len(v) 

    return finalalleles

  def scrape_blast(self,type=""):
    hypo = list()
    type2db = type.capitalize() + 's'
    if type == 'expec':
      type2db = 'Expacs'
    folder = type
    folder = folder.replace('seq_type', 'mlst')
    folder = folder.replace('core_', 'cg')
 
    q_list = glob.glob("{}/blast_search/{}/*".format(self.sampledir, folder))
    organism = self.lims_fetcher.get_organism_refname(self.name)
    if not organism:
      organism = self.lims_fetcher.data['organism']
    self.db_pusher.upd_rec({'CG_ID_sample' : self.name}, 'Samples', {'organism': organism})
    res_cols = self.db_pusher.get_columns('{}'.format(type2db))

    try:
      old_ref = ""
      for file in q_list:
        filename = os.path.basename(file).rsplit('.',1)[0] #Removes suffix
        if filename == 'lactam':
          filename = 'beta-lactam'
        if type == 'resistance':
          ref_file = "{}/{}.fsa".format(self.config["folders"]["resistances"], filename)
        elif type == 'expec':
          ref_file = self.config["folders"]["expec"]
        elif type == 'core_seq_type':
          ref_file = "{}/{}/main.fsa".format(self.config['folders']['cgmlst'], organism)
        elif type == 'seq_type':
          ref_file =  "{}/{}/{}.tfa".format(self.config['folders']['references'], organism, filename.split('_')[-1])

        #Removes a lot of cgmlst inits
        if old_ref != ref_file:
          old_ref = ref_file
          locilengths = self.get_locilengths(ref_file)
          

        with open("{}".format(file), 'r') as sample:
          for line in sample:
            #Ignore commented fields
            if not line[0] == '#':

              elem_list = line.rstrip().split("\t")
              if not elem_list[1] == 'N/A':
                hypo.append(dict())
                hypo[-1]["CG_ID_sample"] = self.name
                hypo[-1]["identity"] = elem_list[4]
                hypo[-1]["evalue"] = elem_list[5]
                hypo[-1]["bitscore"] = elem_list[6]
                if int(elem_list[7]) < int(elem_list[8]):
                  hypo[-1]["contig_start"] = int(elem_list[7])
                  hypo[-1]["contig_end"] = int(elem_list[8])
                else:
                  hypo[-1]["contig_start"] = int(elem_list[8])
                  hypo[-1]["contig_end"] = int(elem_list[7])
                hypo[-1]["subject_length"] =  int(elem_list[11])

                if type == 'resistance':
                  hypo[-1]["instance"] = filename
                  partials = re.search(r'(.+)_(\d+){1,3}(?:_(\w+))*', elem_list[3])
                  hypo[-1]["reference"] = partials.group(3)
                  hypo[-1]["gene"] = partials.group(1)
                  if hypo[-1]["gene"] in self.gene2resistance.keys():
                    hypo[-1]["resistance"] = self.gene2resistance[hypo[-1]["gene"]]
                  else:
                    hypo[-1]["{}".format(type)] = hypo[-1]["instance"].capitalize()
                  hypo[-1]["span"] = float(hypo[-1]["subject_length"])/locilengths['>{}'.format(partials.group(0))]

                elif type == 'expec':
                  hypo[-1]["instance"] = filename
                  #Thanks, precompiled list standards
                  if '>' in elem_list[3]:
                    partials = re.search(r'>*(\w+_\w+\.*\w+).+\((\w+)\).+\((\w+)\)_(\w+)_\[.+\]', elem_list[3])
                  else:
                    partials = re.search(r'(\w+)\(gb\|\w+\)_\((\S+)\)_(.+)_\[(\S+)_.+\]_\[\S+\]', elem_list[3])
                  if not partials:
                    partials = re.search(r'(\w+\.*\w+)\:*\w*_*(?:\(\w+\-\w+\))*_\((\w+)\)_([^[]+)\[\S+\]', elem_list[3])
                  #NC/Protein reference
                  hypo[-1]["reference"] = partials.group(1)
                  #Full gene name
                  hypo[-1]["gene"] = partials.group(2)
                  #More generic group
                  hypo[-1]["instance"] = partials.group(3).strip('_')
                  #Description
                  if len(partials.groups()) >= 4:
                    hypo[-1]["virulence"] = partials.group(4).replace('_', ' ').capitalize()
                  else:
                    hypo[-1]["virulence"] = ""
                  hypo[-1]["span"] = float(hypo[-1]["subject_length"])/locilengths['>{}'.format(elem_list[0])]

                elif type == 'seq_type':
                  partials = re.search(r'(.+)_(\d+){1,3}(?:_(\w+))*', elem_list[3])
                  hypo[-1]["loci"] = partials.group(1)
                  hypo[-1]["allele"] = int(partials.group(2))
                  hypo[-1]["span"] = float(hypo[-1]["subject_length"])/locilengths['>{}'.format(partials.group(0))]

                elif type == 'core_seq_type':
                  partials = re.search(r'(.+)_(\d+){1,3}(?:_(\w+))*', elem_list[2])
                  hypo[-1]["allele"] = int(elem_list[0])
                  hypo[-1]["loci"] = filename
                  hypo[-1]["span"] = float(hypo[-1]["subject_length"])/locilengths['>{}_{}'.format(hypo[-1]["loci"], hypo[-1]["allele"])] 

                # split elem 2 into contig node_NO, length, cov
                nodeinfo = elem_list[2].split('_')
                hypo[-1]["contig_name"] = "{}_{}".format(nodeinfo[0], nodeinfo[1])
                hypo[-1]["contig_length"] = int(nodeinfo[3])
                hypo[-1]["contig_coverage"] = nodeinfo[5]

      self.logger.info("{} candidate {} hits found".format(len(hypo), type2db))
    except Exception as e:
      self.logger.error("Unable to process the pattern of {}".format(str(e)))

    #Cleanup of overlapping hits
    if type == 'seq_type' or type == 'core_seq_type':
      identifier = 'loci'
    elif type == 'resistance' or type == 'expec':
      identifier = 'gene'
    ind = 0
    while ind < len(hypo)-1:
      targ = ind+1
      while targ < len(hypo):
        ignore = False
        if hypo[ind]["contig_name"] == hypo[targ]["contig_name"]:
          #Overlapping or shared gene 
          if (hypo[ind]["contig_start"] >= hypo[targ]["contig_start"] and hypo[ind]["contig_start"] <= hypo[targ]["contig_end"]) or\
              (hypo[ind]["contig_end"] >= hypo[targ]["contig_start"] and hypo[ind]["contig_end"] <= hypo[targ]["contig_end"]) or \
              (hypo[ind][identifier] == hypo[targ][identifier]):
            #Rightmost is worse
            if float(hypo[ind]["identity"])*(1-abs(1-hypo[ind]["span"])) > float(hypo[targ]["identity"])*(1-abs(1-hypo[targ]["span"])):
              del hypo[targ]
              ignore = True
            #Leftmost is worse
            elif float(hypo[ind]["identity"])*(1-abs(1-hypo[ind]['span'])) < float(hypo[targ]["identity"])*(1-abs(1-hypo[targ]['span'])):
              del hypo[ind]
              targ = ind +1
              ignore = True
            #Identical identity and span, seperating based on contig coverage
            else:
              #Rightmost is worse
              if float(hypo[ind]["contig_coverage"]) >= float(hypo[targ]["contig_coverage"]):
                del hypo[targ]
                ignore = True
              #Leftmost is worse
              elif float(hypo[ind]["contig_coverage"]) < float(hypo[targ]["contig_coverage"]):
                del hypo[ind]
                targ = ind +1
                ignore = True
        if not ignore:
          targ += 1
        else:
          pass
      ind += 1

    self.logger.info("{} {} hits were added after removing overlaps and duplicate hits".format( len(hypo), type))
    for hit in hypo:
      #self.logger.info("Kept {}:{} with span {} and id {}".format(hit['loci'],hit["allele"], hit['span'],hit['identity']))
      self.db_pusher.add_rec(hit, '{}'.format(type2db))
  
    if type == 'seq_type':
      try:
        ST = self.db_pusher.alleles2st(self.name)
        self.db_pusher.upd_rec({'CG_ID_sample':self.name}, 'Samples', {'ST':ST})
        self.logger.info("Sample {} received ST {}".format(self.name, ST))
      except Exception as e:
        self.logger.warning("Unable to type sample {} due to data value '{}'".format(self.name, str(e)))

  def load_resistances(self):
    """Legacy function, loads common resistance names for genes from notes file"""
    conversions = dict()
    with open("{}/notes.txt".format(self.config["folders"]["resistances"])) as fh:
      for line in fh:
        if '#' not in line:
          line = line.split(':')
          cropped = re.sub(' resistance', '', line[1])
          conversions[line[0]] = cropped
          #Workaround for case issues
          conversions[line[0].lower()] = cropped
    return conversions

  def scrape_alignment(self):
    """Scrapes a single alignment result"""
    ins_list = list()
    cov_dict = dict()
    align_dict = dict()
    align_dict["reference_genome"] = self.lims_fetcher.data['reference']

    #Reading
    q_list = glob.glob("{}/alignment/*.stats.*".format(self.sampledir))
    map_rate = 0.0
    median_ins = 0
    ref_len = 0.0
    tot_reads = 0
    tot_map = 0
    duprate = 0.0
    for file in q_list:
      with open(file, 'r') as fh:
       type = file.split('.')[-1]
       for line in fh.readlines():
         lsplit = line.rstrip().split('\t')
         if type == 'raw':
           try:
             tot_reads = int(lsplit[0])
           except Exception as e:
             pass
         elif type == 'ins':
           if len(lsplit) >= 18 and lsplit[-12] in ['FF','FR']:
             try:
               median_ins = int(lsplit[0])
             except Exception as e:
               pass
         elif type == 'cov':
           cov_dict[lsplit[1]] = int(lsplit[2])
         elif type == 'ref':
           if lsplit[0] != '*' and len(lsplit) >= 2:
             ref_len = ref_len + int(lsplit[1])
         elif type == 'dup':
           if lsplit[0] == 'Unknown Library':
             try:
               duprate = float(lsplit[8]) 
             except Exception as e:
               duprate = -1.0
         elif type == 'map':
           dsplit = line.rstrip().split(' ')
           if len(dsplit)>= 5 and dsplit[4] == 'total':
             tot_map = int(dsplit[0])
           elif len(dsplit)>=4 and dsplit[3] == 'mapped':
             if tot_map > 0:
               map_rate = int(dsplit[0])/float(tot_map)

    #Mangling
    sumz, plus10, plus30, plus50, plus100, total = 0, 0, 0, 0, 0, 0
    for k, v in cov_dict.items():
      sumz += int(k)*v
      total += v
      if int(k) > 10:
        plus10 += v
      if int(k) > 30:
        plus30 += v
      if int(k) > 50:
        plus50 += v
      if int(k) > 100:
        plus100 += v
    if total > 0:
      align_dict['coverage_10x'] = plus10/float(ref_len)
      align_dict['coverage_30x'] = plus30/float(ref_len)
      align_dict['coverage_50x'] = plus50/float(ref_len)
      align_dict['coverage_100x'] = plus100/float(ref_len)
    else:
      align_dict['coverage_10x'] = 0.0
      align_dict['coverage_30x'] = 0.0
      align_dict['coverage_50x'] = 0.0
      align_dict['coverage_100x'] = 0.0
 
    align_dict['mapped_rate'] = map_rate
    align_dict['insert_size'] = median_ins
    if ref_len > 0:
      align_dict['duplication_rate'] = duprate
      align_dict['average_coverage'] = sumz/float(ref_len)
    else:
      align_dict['duplication_rate'] = 0.0
      align_dict['average_coverage'] = 0.0
    align_dict['total_reads'] = tot_reads
    self.db_pusher.upd_rec({'CG_ID_sample' : self.name}, 'Samples', align_dict)
