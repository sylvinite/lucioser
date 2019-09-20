"""This is the main entry point of microSALT.
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import logging
import json
import os
import re
import subprocess
import sys
import yaml

from pkg_resources import iter_entry_points
from microSALT import __version__, config, wd
from microSALT.utils.scraper import Scraper
from microSALT.utils.job_creator import Job_Creator
from microSALT.utils.reporter import Reporter
from microSALT.utils.referencer import Referencer
from microSALT.store.lims_fetcher import LIMS_Fetcher


if config == '':
  click.echo("ERROR: No properly set-up config under neither envvar MICROSALT_CONFIG nor ~/.microSALT/config.json. Exiting.")
  sys.exit(-1)
else:
  #Makes sure DB inherits correct permissions if freshly created
  bash_cmd="touch {}".format(re.search('sqlite\:\/\/\/(.+)', config['database']['SQLALCHEMY_DATABASE_URI']).group(1))
  proc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
  output, error = proc.communicate()

def done():
  click.echo("Execution finished!")

@click.group()
@click.version_option(__version__)
@click.pass_context
def root(ctx):
  """microbial Sequence Analysis and Loci-based Typing (microSALT) pipeline """
  ctx.obj = {}
  ctx.obj['config'] = config
  logger = logging.getLogger('main_logger')
  logger.setLevel(logging.DEBUG)
  fh = logging.FileHandler(os.path.expanduser(config['folders']['log_file']))
  fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
  logger.addHandler(fh)
  ch = logging.StreamHandler()
  ch.setLevel(logging.INFO)
  ch.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))
  logger.addHandler(ch)
  ctx.obj['log'] = logger
  if os.path.basename(os.path.expandvars('$CONDA_PREFIX'))[0] == 'D':
    ctx.obj['config']['folders']['expec'] = os.path.abspath(os.path.join(wd , '../unique_references/ExPEC.fsa'))
  else:
    ctx.obj['config']['folders']['expec'] = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'expec/ExPEC.fsa'))
  ctx.obj['config']['folders']['adapters'] = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'share/trimmomatic-0.39-1/adapters/'))
  ctx.obj['config']['folders']['st_profiles'] = "{}/ST_profiles".format(ctx.obj['config']['folders']['references'])
  ctx.obj['config']['folders']['st_loci'] = "{}/ST_loci".format(ctx.obj['config']['folders']['references'])
  ctx.obj['config']['folders']['resistances'] = "{}/resfinder_db".format(ctx.obj['config']['folders']['references'])
  ctx.obj['config']['folders']['virulences'] = "{}/virulencefinder_db".format(ctx.obj['config']['folders']['references'])
  ctx.obj['config']['folders']['plasmids'] = "{}/plasmidfinder_db".format(ctx.obj['config']['folders']['references'])
  ctx.obj['config']['folders']['fimhs'] = "{}/fimtyper_db".format(ctx.obj['config']['folders']['references'])
  ctx.obj['config']['folders']['genomes'] = "{}/genomes".format(ctx.obj['config']['folders']['references'])
  ctx.obj['config']['folders']['cgmlst'] = "{}/cgmlst".format(ctx.obj['config']['folders']['references'])
  
  scientist=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])
  scientist.check_connection()

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
  ctx.obj['config']['regex']['mail_recipient'] = email
  trimmed = not untrimmed
  careful = not uncareful
  if config != '':
    if os.path.exists(config):
      try:
        with open(os.path.abspath(config), 'r') as conf:
          ctx.obj['config'] = json.load(conf)
        ctx.obj['config']['config_path'] = os.path.abspath(config)
      except Exception as e:
        pass

  ctx.obj['config']['dry'] = dry
  if input != "":
    project_dir = os.path.abspath(input)
    if not project_id in project_dir:
      click.echo("Path does not contain project id. Exiting.")
      sys.exit(-1)
  else:
    project_dir = "{}/{}".format(ctx.obj['config']['folders']['seqdata'], project_id)
    if not os.path.isdir(project_dir):
      click.echo("Sequence data folder for {} does not exist.".format(project_id))
      sys.exit(-1)

  click.echo("Checking versions of references..")
  try:
    if not skip_update:
      fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
      fixer.identify_new(project_id,project=True)
      fixer.update_refs()
      click.echo("Version check done. Creating sbatch jobs")
    else:
      click.echo("Skipping version check.")
  except Exception as e:
    click.echo("{}".format(e))

  manager = Job_Creator(project_dir, ctx.obj['config'], ctx.obj['log'],trim=trimmed,qc_only=qc_only, careful=careful)
  manager.project_job()
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
  ctx.obj['config']['regex']['mail_recipient'] = email
  trimmed = not untrimmed
  careful = not uncareful
  if config != '':
    if os.path.exists(config):
      try:
        with open(os.path.abspath(config), 'r') as conf:
          ctx.obj['config'] = json.load(conf)
        ctx.obj['config']['config_path'] = os.path.abspath(config)
      except Exception as e:
        pass

  ctx.obj['config']['dry'] = dry
  scientist=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])
  try:
    scientist.load_lims_sample_info(sample_id)
  except Exception as e:
    click.echo("Unable to load LIMS sample info.")

  if input != "":
    sample_dir = os.path.abspath(input)
    if not sample_id in sample_dir:
      click.echo("Path does not contain sample id. Exiting.")
      sys.exit(-1)
  else:
    sample_dir = "{}/{}/{}".format(ctx.obj['config']['folders']['seqdata'], scientist.data['CG_ID_project'] ,sample_id)
    if not os.path.isdir(sample_dir):
      click.echo("Sequence data folder for {} does not exist.".format(sample_id))
      sys.exit(-1)

  click.echo("Checking versions of references..")
  try:
    if not skip_update:
      fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
      fixer.identify_new(sample_id,project=False) 
      fixer.update_refs()
      click.echo("Version check done. Creating sbatch job")
    else:
      click.echo("Skipping version check.")
    worker = Job_Creator(sample_dir, ctx.obj['config'], ctx.obj['log'], trim=trimmed,qc_only=qc_only, careful=careful)
    worker.project_job(single_sample=True)
  except Exception as e:
    click.echo("Unable to process sample {} due to '{}'".format(sample_id,e))
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
  ctx.obj['config']['regex']['mail_recipient'] = email
  trimmed = not untrimmed
  careful = not uncareful
  if config != '':
    if os.path.exists(config):
      try:
        with open(os.path.abspath(config), 'r') as conf:
          ctx.obj['config'] = json.load(conf)
        ctx.obj['config']['config_path'] = os.path.abspath(config)
      except Exception as e:
        pass

  ctx.obj['config']['dry'] = dry
  scientist=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])

  pool = []
  if input != "":
    collection_dir = os.path.abspath(input)
  else:
    collection_dir = "{}/{}".format(ctx.obj['config']['folders']['seqdata'], collection_id)
  for sample in os.listdir(collection_dir):
    if os.path.isdir("{}/{}".format(collection_dir,sample)):
      pool.append(sample)
  if pool == []:
    click.echo("Input collection lacks any valid of samples")
    sys.exit(-1)

  pool_cg = []
  for sample in pool:
    try:
      scientist.load_lims_sample_info(sample,warnings=True)
      pool_cg.append(scientist.data['CG_ID_sample'])
    except Exception as e:
      click.echo("Unable to load LIMS sample info for sample {}.".format(sample))

  click.echo("Checking versions of references..")
  for sample in pool_cg:
    try:
      if not skip_update:
        fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
        fixer.identify_new(sample,project=False)
    except Exception as e:
      click.echo("Unable to update references for sample {} due to '{}'".format(sample,str(e)))
  if not skip_update:
    fixer.update_refs()
    click.echo("Version check done. Creating sbatch job")
  else:
    click.echo("Skipping version check.")

  try: 
    worker = Job_Creator(collection_dir, ctx.obj['config'], ctx.obj['log'], trim=trimmed,qc_only=qc_only, careful=careful, pool=pool_cg)
    worker.project_job()
  except Exception as e:
    click.echo("Unable to process collection due to '{}'".format(str(e)))
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
  if config != '':
    if os.path.exists(config):
      try:
        with open(os.path.abspath(config), 'r') as conf:
          ctx.obj['config'] = json.load(conf)
        ctx.obj['config']['config_path'] = os.path.abspath(config)
      except Exception as e:
        pass

  ctx.obj['config']['rerun'] = rerun
  ctx.obj['config']['regex']['mail_recipient'] = email

  if input != "":
    sample_dir = os.path.abspath(input)
    if not sample_id in sample_dir:
      click.echo("Path does not contain sample id. Exiting.")
      sys.exit(-1)
  else:
    prohits = [x for x in os.listdir(ctx.obj['config']['folders']['results']) if x.startswith("{}_".format(sample_id))]
    if len(prohits) > 1:
      click.echo("Multiple instances of that analysis exists. Specify full path using --input")
      sys.exit(-1)
    elif len(prohits) <1:
      click.echo("No analysis folder prefixed by {} found.".format(sample_id))
      sys.exit(-1)
    else:
      sample_dir = "{}/{}".format(ctx.obj['config']['folders']['results'], prohits[-1])

  scientist=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])
  try:
    scientist.load_lims_sample_info(sample_id)
  except Exception as e:
    click.echo("Unable to load LIMS sample info.")
    sys.exit(-1)

  garbageman = Scraper(sample_dir, ctx.obj['config'], ctx.obj['log'])
  garbageman.scrape_sample()

  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'], scientist.data['CG_ID_project'], output=sample_dir)
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
  if config != '':
    if os.path.exists(config):
      try:
        with open(os.path.abspath(config), 'r') as conf:
          ctx.obj['config'] = json.load(conf)
        ctx.obj['config']['config_path'] = os.path.abspath(config)
      except Exception as e:
        pass

  ctx.obj['config']['rerun'] = rerun
  ctx.obj['config']['regex']['mail_recipient'] = email

  if input != "":
    project_dir = os.path.abspath(input)
    if not project_id in project_dir:
      click.echo("Path does not contain project id. Exiting.")
      sys.exit(-1)
  else:
    prohits = [x for x in os.listdir(ctx.obj['config']['folders']['results']) if x.startswith("{}_".format(project_id))]
    if len(prohits) > 1:
      click.echo("Multiple instances of that analysis exists. Specify full path using --input")
      sys.exit(-1)
    elif len(prohits) <1:
      click.echo("No analysis folder prefixed by {} found.".format(project_id))
      sys.exit(-1)
    else:
      project_dir = "{}/{}".format(ctx.obj['config']['folders']['results'], prohits[-1])

  garbageman = Scraper(project_dir, ctx.obj['config'], ctx.obj['log'])
  garbageman.scrape_project()
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
  if config != '':
    if os.path.exists(config):
      try:
        with open(os.path.abspath(config), 'r') as conf:
          ctx.obj['config'] = json.load(conf)
        ctx.obj['config']['config_path'] = os.path.abspath(config)
      except Exception as e:
        pass

  ctx.obj['config']['rerun'] = False
  ctx.obj['config']['regex']['mail_recipient'] = email

  pool = []
  if input != "":
    collection_dir = os.path.abspath(input)
  else:
    prohits = [x for x in os.listdir(ctx.obj['config']['folders']['results']) if x.startswith("{}_".format(collection_id))]
    if len(prohits) > 1:
      click.echo("Multiple instances of that analysis exists. Specify full path using --input")
      sys.exit(-1)
    elif len(prohits) <1:
      click.echo("No analysis folder prefixed by {} found.".format(project_id))
      sys.exit(-1)
    else:
      collection_dir = "{}/{}".format(ctx.obj['config']['folders']['results'], prohits[-1])

  for sample in os.listdir(collection_dir):
    if os.path.isdir("{}/{}".format(collection_dir,sample)):
      pool.append(sample)
  if pool == []:
    click.echo("Input collection lacks any valid of samples")
    sys.exit(-1)

  pool_cg = []
  scientist=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])
  for sample in pool:
    try:
      scientist.load_lims_sample_info(sample,warnings=True)
      pool_cg.append(scientist.data['CG_ID_sample'])
    except Exception as e:
      click.echo("Unable to load LIMS sample info for sample {}.".format(sample))
      sys.exit(-1)

  for sample in pool:
    garbageman = Scraper("{}/{}".format(collection_dir, sample), ctx.obj['config'], ctx.obj['log'])
    garbageman.scrape_sample()
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
    sys.exit(-1)
  click.echo("Checking versions of all references..")
  referee = Referencer(ctx.obj['config'], ctx.obj['log'],force=force)
  referee.update_refs()

@refer.command()
@click.pass_context
def list(ctx):
  """ Lists all stored organisms """
  refe = Referencer(ctx.obj['config'], ctx.obj['log'])
  click.echo("Currently stored organisms:")
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

@resync.command()
@click.option('--type', default='html', type=click.Choice(['report', 'list']), help="Output format")
@click.option('--customer', default='all', help="Customer id filter")
@click.option('--skip_update', default=False, help="Skips downloading of references", is_flag=True)
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.pass_context
def review(ctx, type, customer, skip_update, email):
  """Generates information about novel ST"""
  #Trace exists by some samples having pubMLST_ST filled in. Make trace function later
  ctx.obj['config']['regex']['mail_recipient'] = email
  fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
  if not skip_update:
    fixer.update_refs()
    fixer.resync()
  click.echo("Version check done. Generating output")
  if type=='report':
    codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'])
    codemonkey.report(type='st_update', customer=customer)
  elif type=='list':
    fixer.resync(type=type)
  done()

@resync.command()
@click.argument('sample_name')
@click.pass_context
def overwrite(ctx,sample_name):
  """Flags sample as resolved"""
  fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
  fixer.resync(type='overwrite', sample=sample_name)
  done()
