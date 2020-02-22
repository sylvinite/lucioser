"""This is the main entry point of microSALT.
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import json
import os
import re
import subprocess
import sys
import yaml

from pkg_resources import iter_entry_points
from microSALT import __version__, config, logger, wd
from microSALT.utils.scraper import Scraper
from microSALT.utils.job_creator import Job_Creator
from microSALT.utils.reporter import Reporter
from microSALT.utils.referencer import Referencer
from microSALT.store.lims_fetcher import LIMS_Fetcher

if config == '':
  click.echo("ERROR - No properly set-up config under neither envvar MICROSALT_CONFIG nor ~/.microSALT/config.json. Exiting.")
  click.abort()

def set_cli_config(config):
  if config != '':
    if os.path.exists(config):
      try:
        t = ctx.obj['config']
        with open(os.path.abspath(config), 'r') as conf:
          ctx.obj['config'] = json.load(conf)
        ctx.obj['config']['folders']['expec'] = t['folders']['expec']
        ctx.obj['config']['folders']['adapters'] = t['folders']['adapters']
        ctx.obj['config']['config_path'] = os.path.abspath(config)
      except Exception as e:
        pass

def done():
  click.echo("INFO - Execution finished!")

@click.group()
@click.version_option(__version__)
@click.pass_context
def root(ctx):
  """microbial Sequence Analysis and Loci-based Typing (microSALT) pipeline """
  ctx.obj = {}
  ctx.obj['config'] = config
  ctx.obj['log'] = logger
  lims_obj=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])
  lims_obj.check_connection() 

@root.group()
@click.pass_context
def analyse(ctx):
  """Basic sequence and resistance typing"""
  pass

@root.group()
@click.pass_context
def utils(ctx):
  """ Utilities for specific purposes """
  pass

@utils.group()
@click.pass_context
def finish(ctx):
  """Reupload typing analysis and generate results"""
  pass

@utils.group()
@click.pass_context
def refer(ctx):
  """ Manipulates MLST organisms """
  pass

@analyse.command()
@click.argument('project_id')
@click.option('--input', help='Full path to project folder',default="")
@click.option('--dry', help="Builds instance without posting to SLURM", default=False, is_flag=True)
@click.option('--qc_only', help="Only runs QC (alignment stats)", default=False, is_flag=True)
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--skip_update', default=False, help="Skips downloading of references", is_flag=True)
@click.option('--untrimmed', help="Use untrimmed input data", default=False, is_flag=True)
@click.option('--uncareful', help="Avoids running SPAdes in careful mode. Sometimes fix assemblies", default=False, is_flag=True)
@click.pass_context
def project(ctx, project_id, input, dry, config, email, qc_only, untrimmed, skip_update, uncareful):
  """Analyse a whole project"""
  trimmed = not untrimmed
  careful = not uncareful
  set_cli_config(config)
  ctx.obj['config']['regex']['mail_recipient'] = email
  ctx.obj['config']['dry'] = dry

  if input != "":
    project_dir = os.path.abspath(input)
    if not project_id in project_dir:
      click.echo("ERROR - Path does not contain project id. Exiting.")
      click.abort()
  else:
    project_dir = "{}/{}".format(ctx.obj['config']['folders']['seqdata'], project_id)
    if not os.path.isdir(project_dir):
      click.echo("ERROR - Sequence data folder for {} does not exist.".format(project_id))
      click.abort()

  ext_refs = Referencer(ctx.obj['config'], ctx.obj['log'])
  click.echo("INFO - Checking versions of references..")
  try:
    if not skip_update:
      ext_refs.identify_new(project_id,project=True)
      ext_refs.update_refs()
      click.echo("INFO - Version check done. Creating sbatch jobs")
    else:
      click.echo("INFO - Skipping version check.")
  except Exception as e:
    click.echo("{}".format(e))

  run_creator = Job_Creator(project_dir, ctx.obj['config'], ctx.obj['log'],trim=trimmed,qc_only=qc_only, careful=careful)
  run_creator.project_job()
  done() 

@analyse.command()
@click.argument('sample_id')
@click.option('--input', help='Full path to sample folder', default="")
@click.option('--dry', help="Builds instance without posting to SLURM", default=False, is_flag=True)
@click.option('--qc_only', help="Only runs QC (alignment stats)", default=False, is_flag=True)
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--untrimmed', help="Use untrimmed input data", default=False, is_flag=True)
@click.option('--skip_update', default=False, help="Skips downloading of references", is_flag=True)
@click.option('--uncareful', help="Avoids running SPAdes in careful mode. Sometimes fix assemblies", default=False, is_flag=True)
@click.pass_context
def sample(ctx, sample_id, input, dry, config, email, qc_only, untrimmed, skip_update, uncareful):
  """Analyse a single sample"""

  trimmed = not untrimmed
  careful = not uncareful
  set_cli_config(config)
  ctx.obj['config']['regex']['mail_recipient'] = email
  ctx.obj['config']['dry'] = dry

  lims_obj=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])
  try:
    lims_obj.load_lims_sample_info(sample_id)
  except Exception as e:
    click.echo("ERROR - Unable to load LIMS sample info.")

  if input != "":
    sample_dir = os.path.abspath(input)
    if not sample_id in sample_dir:
      click.echo("ERROR - Path does not contain sample id. Exiting.")
      click.abort()
  else:
    sample_dir = "{}/{}/{}".format(ctx.obj['config']['folders']['seqdata'], lims_obj.data['CG_ID_project'] ,sample_id)
    if not os.path.isdir(sample_dir):
      click.echo("ERROR - Sequence data folder for {} does not exist.".format(sample_id))
      click.abort()

  ext_refs = Referencer(ctx.obj['config'], ctx.obj['log'])
  click.echo("INFO - Checking versions of references..")
  try:
    if not skip_update:
      ext_refs.identify_new(sample_id,project=False) 
      ext_refs.update_refs()
      click.echo("INFO - Version check done. Creating sbatch job")
    else:
      click.echo("INFO - Skipping version check.")
    worker = Job_Creator(sample_dir, ctx.obj['config'], ctx.obj['log'], trim=trimmed,qc_only=qc_only, careful=careful)
    worker.project_job(single_sample=True)
  except Exception as e:
    click.echo("ERROR - Unable to process sample {} due to '{}'".format(sample_id,e))
  done()


@analyse.command()
@click.argument('collection_id')
@click.option('--input', help='Full path to sample folder', default="")
@click.option('--dry', help="Builds instance without posting to SLURM", default=False, is_flag=True)
@click.option('--qc_only', help="Only runs QC (alignment stats)", default=False, is_flag=True)
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--untrimmed', help="Use untrimmed input data", default=False, is_flag=True)
@click.option('--skip_update', default=False, help="Skips downloading of references", is_flag=True)
@click.option('--uncareful', help="Avoids running SPAdes in careful mode. Sometimes fix assemblies", default=False, is_flag=True)
@click.pass_context
def collection(ctx, collection_id, input, dry, qc_only, config, email, untrimmed, skip_update, uncareful):
  """Analyse a collection of samples"""
  trimmed = not untrimmed
  careful = not uncareful
  set_cli_config(config)
  ctx.obj['config']['regex']['mail_recipient'] = email
  ctx.obj['config']['dry'] = dry

  lims_obj=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])

  pool = []
  if input != "":
    collection_dir = os.path.abspath(input)
    if not collection_id in collection_dir:
      click.echo("ERROR - Path does not contain collection id. Exiting.")
      click.abort()
  else:
    collection_dir = "{}/{}".format(ctx.obj['config']['folders']['seqdata'], collection_id)
    if not os.path.isdir(collection_dir):
      click.echo("Collection data folder does not exist. Exiting.")
      click.abort()

  for sample in os.listdir(collection_dir):
    if os.path.isdir("{}/{}".format(collection_dir,sample)):
      pool.append(sample)
  #if pool == []:
  #  click.echo("Input collection lacks any valid of samples")
  #  click.abort()

  pool_cg = []
  for sample in pool:
    try:
      lims_obj.load_lims_sample_info(sample,warnings=True)
      pool_cg.append(lims_obj.data['CG_ID_sample'])
    except Exception as e:
      click.echo("ERROR - Unable to load LIMS sample info for sample {}.".format(sample))

  ext_refs = Referencer(ctx.obj['config'], ctx.obj['log'])
  click.echo("INFO - Checking versions of references..")
  for sample in pool_cg:
    try:
      if not skip_update:
        ext_refs.identify_new(sample,project=False)
    except Exception as e:
      click.echo("ERROR - Unable to update references for sample {} due to '{}'".format(sample,str(e)))
  if not skip_update:
    ext_refs.update_refs()
    click.echo("INFO - Version check done. Creating sbatch job")
  else:
    click.echo("INFO - Skipping version check")

  try: 
    worker = Job_Creator(collection_dir, ctx.obj['config'], ctx.obj['log'], trim=trimmed,qc_only=qc_only, careful=careful, pool=pool_cg)
    worker.project_job()
  except Exception as e:
    click.echo("ERROR - Unable to process collection due to '{}'".format(str(e)))
  done()

@finish.command()
@click.argument('sample_id')
@click.option('--rerun', is_flag=True, default=False, help='Overwrite existing data')
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--input', help='Full path to result sample folder', default="")
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--report', default='default', type=click.Choice(['default', 'qc']))
@click.pass_context
def sample(ctx, sample_id, rerun, email, input, config, report):
  """Parse results from analysing a single sample"""
  set_cli_config(config)
  ctx.obj['config']['regex']['mail_recipient'] = email
  ctx.obj['config']['rerun'] = rerun

  if input != "":
    sample_dir = os.path.abspath(input)
    if not sample_id in sample_dir:
      click.echo("ERROR - Path does not contain sample id. Exiting.")
      click.abort()
  else:
    prohits = [x for x in os.listdir(ctx.obj['config']['folders']['results']) if x.startswith("{}_".format(sample_id))]
    if len(prohits) > 1:
      click.echo("ERROR - Multiple instances of that analysis exists. Specify full path using --input")
      click.abort()
    elif len(prohits) <1:
      click.echo("ERROR - No analysis folder prefixed by {} found.".format(sample_id))
      click.abort()
    else:
      sample_dir = "{}/{}".format(ctx.obj['config']['folders']['results'], prohits[-1])

  lims_obj=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])
  try:
    lims_obj.load_lims_sample_info(sample_id)
  except Exception as e:
    click.echo("ERROR - Unable to load LIMS sample info.")
    click.abort()

  res_scraper = Scraper(sample_dir, ctx.obj['config'], ctx.obj['log'])
  res_scraper.scrape_sample()

  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'], lims_obj.data['CG_ID_project'], output=sample_dir)
  codemonkey.report(report)
  done()

@finish.command()
@click.argument('project_id')
@click.option('--rerun', is_flag=True, default=False, help='Overwrite existing data')
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--input', help='Full path to result project folder', default="")
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--report', default='default', type=click.Choice(['default', 'qc']))
@click.pass_context
def project(ctx, project_id, rerun, email, input, config, report):
  """Parse results from analysing a single project"""
  set_cli_config(config)
  ctx.obj['config']['regex']['mail_recipient'] = email
  ctx.obj['config']['rerun'] = rerun

  if input != "":
    project_dir = os.path.abspath(input)
    if not project_id in project_dir:
      click.echo("ERROR - Path does not contain project id. Exiting.")
      click.abort()
  else:
    prohits = [x for x in os.listdir(ctx.obj['config']['folders']['results']) if x.startswith("{}_".format(project_id))]
    if len(prohits) > 1:
      click.echo("ERROR - Multiple instances of that analysis exists. Specify full path using --input")
      click.abort()
    elif len(prohits) <1:
      click.echo("ERROR - No analysis folder prefixed by {} found.".format(project_id))
      click.abort()
    else:
      project_dir = "{}/{}".format(ctx.obj['config']['folders']['results'], prohits[-1])

  res_scraper = Scraper(project_dir, ctx.obj['config'], ctx.obj['log'])
  res_scraper.scrape_project()
  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'], project_id, output=project_dir)
  codemonkey.report(report)
  done()

@finish.command()
@click.argument('collection_id')
@click.option('--rerun', is_flag=True, default=False, help='Overwrite existing data')
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--input', help='Full path to result sample folder', default="")
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--report', default='default', type=click.Choice(['default', 'qc','motif_overview']))
@click.pass_context
def collection(ctx, collection_id, rerun, email, input, config, report):
  """Parse results from analysing a set of sample"""
  set_cli_config(config)
  ctx.obj['config']['regex']['mail_recipient'] = email
  ctx.obj['config']['rerun'] = rerun

  pool = []
  if input != "":
    collection_dir = os.path.abspath(input)
  else:
    prohits = [x for x in os.listdir(ctx.obj['config']['folders']['results']) if x.startswith("{}_".format(collection_id))]
    if len(prohits) > 1:
      click.echo("ERROR - Multiple instances of that analysis exists. Specify full path using --input")
      click.abort()
    elif len(prohits) <1:
      click.echo("ERROR - No analysis folder prefixed by {} found.".format(project_id))
      click.abort()
    else:
      collection_dir = "{}/{}".format(ctx.obj['config']['folders']['results'], prohits[-1])

  for sample in os.listdir(collection_dir):
    if os.path.isdir("{}/{}".format(collection_dir,sample)):
      pool.append(sample)
  if pool == []:
    click.echo("ERROR - Input collection lacks any valid of samples")
    click.abort()

  pool_cg = []
  lims_obj=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])
  for sample in pool:
    try:
      lims_obj.load_lims_sample_info(sample,warnings=True)
      pool_cg.append(lims_obj.data['CG_ID_sample'])
    except Exception as e:
      click.echo("ERROR - Unable to load LIMS sample info for sample {}.".format(sample))
      click.abort()

  for sample in pool:
    res_scraper = Scraper("{}/{}".format(collection_dir, sample), ctx.obj['config'], ctx.obj['log'])
    res_scraper.scrape_sample()
  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'], collection_id, output=collection_dir, collection=True)
  codemonkey.report(report)
  done()


@refer.command()
@click.argument('organism')
@click.pass_context
@click.option('--force', help="Redownloads existing organism", default=False, is_flag=True)
def add(ctx, organism, force):
  """ Adds a new internal organism from pubMLST """
  referee = Referencer(ctx.obj['config'], ctx.obj['log'],force=force)
  try:
    referee.add_pubmlst(organism)
  except Exception as e:
    click.echo(e.args[0])
    click.abort()
  click.echo("INFO - Checking versions of all references..")
  referee = Referencer(ctx.obj['config'], ctx.obj['log'],force=force)
  referee.update_refs()

@refer.command()
@click.pass_context
def list(ctx):
  """ Lists all stored organisms """
  refe = Referencer(ctx.obj['config'], ctx.obj['log'])
  click.echo("INFO - Currently stored organisms:")
  for org in sorted(refe.existing_organisms()):
    click.echo(org.replace("_"," ").capitalize())

@utils.command()
@click.argument('project_name')
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--type', default='default', type=click.Choice(['default', 'typing', 'motif_overview', 'qc', 'json_dump', 'st_update']))
@click.option('--output',help='Full path to output folder',default="")
@click.option('--collection',default=False, is_flag=True)
@click.pass_context
def report(ctx, project_name, email, type, output, collection):
  """Re-generates report for a project"""
  ctx.obj['config']['regex']['mail_recipient'] = email
  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'], project_name, output, collection=collection)
  codemonkey.report(type)
  done()

@utils.command()
@click.pass_context
def view(ctx):
  """Starts an interactive webserver for viewing"""
  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'])
  codemonkey.start_web()

@utils.group()
@click.pass_context
def resync(ctx):
  """Updates internal ST with pubMLST equivalent"""

@utils.command()
@click.option('--dry', help="Builds instance without posting to SLURM", default=False, is_flag=True)
@click.option('--skip_update', default=False, help="Skips downloading of references", is_flag=True)
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.pass_context
def autobatch(ctx, dry, skip_update, email):
  """Analyses all currently unanalysed projects in the seqdata folder"""
  #Trace exists by some samples having pubMLST_ST filled in. Make trace function later
  ctx.obj['config']['regex']['mail_recipient'] = email
  ext_refs = Referencer(ctx.obj['config'], ctx.obj['log'])
  if not skip_update:
    ext_refs.update_refs()
    ext_refs.resync()

  process = subprocess.Popen('squeue --format="%50j" -h -r'.split(), stdout=subprocess.PIPE)
  run_jobs, error = process.communicate()
  run_jobs = run_jobs.splitlines()
  run_jobs = [ re.search("\"(.+)\"", str(x)).group(1).replace(' ', '') for x in run_jobs];
  old_jobs = os.listdir(ctx.obj['config']['folders']['results'])

  for f in os.listdir(ctx.obj['config']['folders']['seqdata']):
      #Project name not found in slurm list 
      if len([s for s in run_jobs if f in s]) == 0:
        #Project name not found in results directories
        if len([s for s in old_jobs if f in s]) == 0:
          if dry:
            click.echo("DRY - microSALT analyse {} --skip_update".format(f))
          else:
            process = subprocess.Popen("microSALT analyse project {} --skip_update".format(f).split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
        elif dry:
          click.echo("INFO - Skipping {} due to existing analysis in results folder".format(f))
      elif dry:
        click.echo("INFO - Skipping {} due to concurrent SLURM run".format(f))

@resync.command()
@click.option('--type', default='list', type=click.Choice(['report', 'list']), help="Output format")
@click.option('--customer', default='all', help="Customer id filter")
@click.option('--skip_update', default=False, help="Skips downloading of references", is_flag=True)
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.pass_context
def review(ctx, type, customer, skip_update, email):
  """Generates information about novel ST"""
  #Trace exists by some samples having pubMLST_ST filled in. Make trace function later
  ctx.obj['config']['regex']['mail_recipient'] = email
  ext_refs = Referencer(ctx.obj['config'], ctx.obj['log'])
  if not skip_update:
    ext_refs.update_refs()
    ext_refs.resync()
  click.echo("INFO - Version check done. Generating output")
  if type=='report':
    codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'])
    codemonkey.report(type='st_update', customer=customer)
  elif type=='list':
    ext_refs.resync(type=type)
  done()

@resync.command()
@click.argument('sample_name')
@click.option('--force', default=False, is_flag=True, help="Resolves sample without checking for pubMLST match")
@click.pass_context
def overwrite(ctx,sample_name, force):
  """Flags sample as resolved"""
  ext_refs = Referencer(ctx.obj['config'], ctx.obj['log'])
  ext_refs.resync(type='overwrite', sample=sample_name, ignore=force)
  done()

@util.command()
@click.argument('sample_id')
@click.option('--input', help='Full path to result sample folder', default="")
@click.option('--dry', help="Builds instance without posting to SLURM", default=False, is_flag=True)
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.pass_context
def contamination(ctx, sample_id, input, dry, config, email):
  """Performs contamination analysis"""
  print("Hello!")
  ctx.obj['config']['regex']['mail_recipient'] = email
  if config != '':
    try:
      with open(os.path.abspath(config), 'r') as conf:
        ctx.obj['config'] = json.load(conf)
    except Exception as e:
      pass

  ctx.obj['config']['dry'] = dry  
  
  ##Scientists stuff here
  
  if input != "":
    sample_dir = os.path.abspath(input)
    if not sample_id in sample_dir:
      print("Path does not contain sample id. Exiting.")
      click.abort()
  else:
    hits = 0
    for i in os.listdir(ctx.obj['config']['folders']['results']):
      if '{}_'.format(sample_id) in i:
        hits = hits+1
        fname = i
    if hits > 1: #Doublechecks only 1 analysis exists
      print("Multiple instances of that analysis exists. Specify full path using --input")
      click.abort()
    elif hits < 1:
      print("No analysis folder prefixed by {} found.".format(sample_id))
      click.abort()
    else:
      sample_dir = "{}/{}".format(ctx.obj['config']['folders']['results'], fname)

  print(sample_dir)

  lims_obj=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])
  
  lims_obj.load_lims_sample_info(sample_id)
  
  #sample_dir = "{}/{}/{}".format(ctx.obj['config']['folders']['seqdata'], lims_obj.data['CG_ID_project'] ,sample_id)
  worker = Job_Creator(sample_dir, ctx.obj['config'], ctx.obj['log'], sample_dir)
  worker.contamination_job() ##supply database etc here  
