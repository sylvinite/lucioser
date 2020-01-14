"""Generates various reports by tapping into Flask and mySQL
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python
import json
import requests
import time
import os
import socket
import sys
import smtplib

from datetime import datetime
from shutil import copyfile

from os.path import basename
from email.mime.text  import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication

from multiprocessing import Process

from microSALT import __version__
from microSALT.server.views import app, session, gen_reportdata, gen_collectiondata
from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.store.orm_models import Samples
from microSALT.store.lims_fetcher import LIMS_Fetcher

class Reporter():

  def __init__(self, config, log, name = "", output = "", collection=False):
    self.db_pusher=DB_Manipulator(config, log)
    self.name = name
    self.collection = collection
    if output == "":
      self.output = os.getcwd()
    else:
      self.output = output
    self.config = config
    self.logger = log
    for k, v in config.items():
      app.config[k] = v
    self.server = Process(target=app.run)
    self.ticketFinder = LIMS_Fetcher(self.config, self.logger)
    self.attachments = list()
    self.filelist = list()
    self.error = False
    self.dt = datetime.now()
    self.filelist = list() 
    self.now = time.strftime("{}.{}.{}_{}.{}.{}".\
    format(self.dt.year, self.dt.month, self.dt.day, self.dt.hour, self.dt.minute, self.dt.second))

  def report(self, type='default', customer='all'):
    self.gen_version(self.name)
    self.start_web()
    if type == 'json_dump':
      self.gen_json()
    elif type == 'default':
      self.gen_typing()
      self.gen_qc()
      self.gen_json(silent=True)
    elif type == 'typing':
      self.gen_typing()
    elif type == 'motif_overview':
      self.gen_motif(motif="resistance")
      self.gen_motif(motif="expec")
    elif type == 'qc':
      self.gen_qc()
    elif type == 'st_update':
      self.gen_STtracker(customer)
    else:
      raise Exception("Report function recieved invalid format")
    self.kill_flask()
    self.mail()
    if self.output == "" or self.output == os.getcwd():
      for file in self.filelist:
        os.remove(file)

  def gen_version(self, name):
    self.db_pusher.get_report(name)
    self.db_pusher.set_report(name)

  def gen_STtracker(self, customer="all", silent=False):
    self.name ="Sequence Type Update"
    try:
      r = requests.get("http://127.0.0.1:5000/microSALT/STtracker/{}".format(customer), allow_redirects=True)
      outname = "{}/ST_updates_{}.html".format(self.output, self.now)
      open(outname, 'wb').write(r.content.decode("iso-8859-1").encode("utf8"))
      self.filelist.append(outname)
      if not silent:
        self.attachments.append(outname)
    except Exception as e:
      self.logger.error("Flask instance currently occupied. Possible rogue process. Retry command")
      self.error = True

  def gen_qc(self,silent=False):
    try:
      self.ticketFinder.load_lims_project_info(self.name)
      last_version = self.db_pusher.get_report(self.name).version
    except Exception as e:
      self.logger.error("Project {} does not exist".format(self.name))
      self.kill_flask()
      sys.exit(-1)
    try:
      q = requests.get("http://127.0.0.1:5000/microSALT/{}/qc".format(self.name), allow_redirects=True)
      outfile = "{}_QC_{}.html".format(self.ticketFinder.data['Customer_ID_project'], last_version)
      output ="{}/{}".format(self.output, outfile)
      storage = "{}/{}".format(self.config['folders']['reports'], outfile)

      if not os.path.isfile(output):
        self.filelist.append(output)
      open(output, 'wb').write(q.content.decode("iso-8859-1").encode("utf8"))
      copyfile(output, storage)

      if not silent:
        self.attachments.append(output)  
    except Exception as e:
      self.logger.error("Flask instance currently occupied. Possible rogue process. Retry command")
      self.error = True
 
  def gen_typing(self,silent=False):
    try:
      self.ticketFinder.load_lims_project_info(self.name)
      last_version = self.db_pusher.get_report(self.name).version
    except Exception as e:
      self.logger.error("Project {} does not exist".format(self.name))
      self.kill_flask()
      sys.exit(-1)
    try:
      r = requests.get("http://127.0.0.1:5000/microSALT/{}/typing/all".format(self.name), allow_redirects=True)
      outfile = "{}_Typing_{}.html".format(self.ticketFinder.data['Customer_ID_project'], last_version)
      output ="{}/{}".format(self.output, outfile)
      storage = "{}/{}".format(self.config['folders']['reports'], outfile)

      if not os.path.isfile(output):
        self.filelist.append(output)
      open(output, 'wb').write(r.content.decode("iso-8859-1").encode("utf8"))
      copyfile(output, storage)
    
      if not silent:
        self.attachments.append(output)
    except Exception as e:
      self.logger.error("Flask instance currently occupied. Possible rogue process. Retry command")
      self.error = True

  def gen_motif(self, motif="resistance", silent=False):
    if motif not in ["resistance", "expec"]:
      self.logger.error("Invalid motif type specified for gen_motif function")
    if self.collection:
      sample_info = gen_collectiondata(self.name)
    else:
      self.ticketFinder.load_lims_project_info(self.name)
      sample_info = gen_reportdata(self.name)
    output = "{}/{}_{}_{}.csv".format(self.output,self.name,motif,self.now)
    if not os.path.isfile(output):
      self.filelist.append(output)
    excel = open(output, "w+")
    motifdict = dict()

    #Load motif & gene names into dict
    for s in sample_info['samples']:
      if motif=='resistance':
        for r in s.resistances:
          if not (r.resistance in motifdict.keys()) and r.threshold == 'Passed':
            motifdict[r.resistance] =list()
          if r.threshold == 'Passed' and not r.gene in motifdict[r.resistance]:
            motifdict[r.resistance].append(r.gene)
      elif motif=='expec':
        for r in s.expacs:
          if not (r.virulence in motifdict.keys()) and r.threshold == 'Passed':
            motifdict[r.virulence] =list()
          if r.threshold == 'Passed' and not r.gene in motifdict[r.virulence]:
            motifdict[r.virulence].append(r.gene)
    for k, v in motifdict.items():
      motifdict[k] = sorted(v)

    #Top 2 Header
    sepfix = "sep=,"
    topline = "Identity {}% & Span {}%,,,".format(self.config['threshold']['motif_id'], self.config['threshold']['motif_span'])
    botline = "CG Sample ID,Sample ID,Organism,Sequence Type,Thresholds"
    for k in sorted(motifdict.keys()):
      genes = [''] * len(motifdict[k])
      active_gene = k.replace(',',' &')
      if active_gene == "":
        active_gene = "Uncategorized hits"
      geneholder = ','.join(genes)
      topline += ",,{}{}".format(active_gene,geneholder)
      resnames = ','.join(sorted(motifdict[k]))
      botline += ",,{}".format(resnames)
    excel.write("{}\n".format(sepfix))
    excel.write("{}\n".format(topline))
    excel.write("{}\n".format(botline))

    #Create each individual row past the 2nd, per iteration
    for s in sample_info['samples']:
      rowdict = dict()
      pref = "{},{},{},{},{}".format(s.CG_ID_sample,s.Customer_ID_sample, s.organism, s.ST_status.replace(',',';'), s.threshold)
      #Load single sample
      if motif=='resistance':
        for r in s.resistances:
          if not (r.resistance in rowdict.keys()) and r.threshold == 'Passed':
            rowdict[r.resistance] =dict()
          if r.threshold == 'Passed' and not r.gene in rowdict[r.resistance]:
            rowdict[r.resistance][r.gene] = r.identity
      elif motif=="expec":
        for r in s.expacs:
          if not (r.virulence in rowdict.keys()) and r.threshold == 'Passed':
            rowdict[r.virulence] =dict()
          if r.threshold == 'Passed' and not r.gene in rowdict[r.virulence]:
            rowdict[r.virulence][r.gene] = r.identity
      #Compare single sample to all
      hits = ""
      for res in sorted(motifdict.keys()):
        if res in rowdict.keys():
          hits += ",1"
          for gen in sorted(motifdict[res]):
            hits += ","
            if gen in rowdict[res].keys():
              #UPD: Change this to identity of hit
              hits +="{}".format(rowdict[res][gen])
            else:
              hits +="0"
        else:
          #Commas eq to res + gen length
          hits += ",0,0"
          pad = ['0'] * len(motifdict[res])
          hits += ','.join(pad)

      excel.write("{}{}\n".format(pref, hits))

    excel.close()
    if not silent:
      self.attachments.append(output)

  def gen_json(self, silent=False):
    report = dict()
    output = "{}/{}_{}.json".format(self.output, self.name, self.now)

    sample_info = gen_reportdata(self.name)
    analyses = ['blast_pubmlst', 'quast_assembly', 'blast_resfinder_resistence', 'picard_markduplicate', 'microsalt_samtools_stats']
    for s in sample_info['samples']:
      report[s.CG_ID_sample] = dict()
      for a in analyses:
        if a == 'blast_resfinder_resistence':
          report[s.CG_ID_sample][a] = list()
        else:
          report[s.CG_ID_sample][a] = dict()

      report[s.CG_ID_sample]['blast_pubmlst'] = {'sequence_type':s.ST_status, 'thresholds':s.threshold}
      report[s.CG_ID_sample]['quast_assembly'] = {'estimated_genome_length':s.genome_length, 'gc_percentage':float(s.gc_percentage), 'n50':s.n50, 'necessary_contigs':s.contigs}
      report[s.CG_ID_sample]['picard_markduplicate'] = {'insert_size':s.insert_size, 'duplication_rate':s.duplication_rate}
      report[s.CG_ID_sample]['microsalt_samtools_stats'] = {'total_reads':s.total_reads, 'mapped_rate':s.mapped_rate,\
                                                            'average_coverage':s.average_coverage, \
                                                            'coverage_10x':s.coverage_10x, 'coverage_30x':s.coverage_30x, \
                                                            'coverage_50x':s.coverage_50x, 'coverage_100x':s.coverage_100x}

      for r in s.resistances:
        if not (r.gene in report[s.CG_ID_sample]['blast_resfinder_resistence']) and r.threshold == 'Passed':
          report[s.CG_ID_sample]['blast_resfinder_resistence'].append(r.gene)

    #json.dumps(report) #Dumps the json directly
    if not os.path.isfile(output):
      self.filelist.append(output)
    with open(output, 'w') as outfile:
      json.dump(report, outfile)
    if not silent:
      self.attachments.append(output)

  def mail(self):
    file_name = self.attachments
    msg = MIMEMultipart()
    if not self.error:
      msg['Subject'] = '{} ({}) Reports'.format(self.name, file_name[0].split('_')[0])
    else:
      msg['Subject'] = '{} Failed Generating Report'.format(self.name)

    if '.' in socket.gethostname():
      msg['From'] = 'microSALT@{}'.format(socket.gethostname())
    else:
      msg['From'] = 'microSALT@{}.com'.format(socket.gethostname()) 
    
    msg['To'] = self.config['regex']['mail_recipient']
   
    if not self.error: 
      for file in self.attachments:
        part = MIMEApplication(open(file).read())  
        part.add_header('Content-Disposition', 'attachment; filename="%s"' % os.path.basename(file))
        msg.attach(part)

    s = smtplib.SMTP('localhost')
    s.connect()
    s.sendmail(msg['From'], msg['To'], msg.as_string())
    s.quit()
    self.logger.info("Mail containing report sent to {} from {}".format(msg['To'], msg['From'])) 

  def start_web(self):
    self.server.start()
    self.logger.info("Started webserver on http://127.0.0.1:5000/")
    #Hinders requests before server goes up
    time.sleep(0.15)

  def kill_flask(self):
    self.server.terminate()
    self.server.join()
    self.logger.info("Closed webserver on http://127.0.0.1:5000/")
