import collections
import logging
import json
import os
import pathlib
import sys

from flask import Flask
from distutils.sysconfig import get_python_lib

__version__ = '2.8.26'

app = Flask(__name__, template_folder='server/templates')
app.config.setdefault('SQLALCHEMY_DATABASE_URI', 'sqlite:///:memory:')
app.config.setdefault('SQLALCHEMY_BINDS', None)
app.config.setdefault('SQLALCHEMY_TRACK_MODIFICATIONS', False)

#Keep track of microSALT installation
wd=os.path.dirname(os.path.realpath(__file__))

# Load configuration
config = ''
logger = ''
default = os.path.join(os.environ['HOME'], '.microSALT/config.json')

if 'MICROSALT_CONFIG' in os.environ:
  try:
    envvar = os.environ['MICROSALT_CONFIG']
    with open(envvar, 'r') as conf:
      config = json.load(conf)
  except Exception as e:
    print("Config error: {}".format(str(e)))
    pass
elif os.path.exists(default):
  try:
    with open(os.path.abspath(default), 'r') as conf:
      config = json.load(conf)
  except Exception as e:
    print("Config error: {}".format(str(e))) 
    pass

# Config dependent section:
if config != '':
  try:
    #Load flask info
    app.config.update(config['database'])

    #Add extrapaths to config
    config['folders']['expec'] = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'unique_references/ExPEC.fsa'))
    #Check if release install exists
    for entry in os.listdir(get_python_lib()):
      if 'microSALT-' in entry:
        config['folders']['expec'] = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'expec/ExPEC.fsa'))
        break
    config['folders']['adapters'] = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'share/trimmomatic-0.39-1/adapters/'))

    #Initialize logger
    logger = logging.getLogger('main_logger')
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(os.path.expanduser(config['folders']['log_file']))
    fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))
    logger.addHandler(ch)

    #Create paths mentioned in config
    for entry in config.keys():
      if entry != '_comment':
        if isinstance(config[entry], str) and '/' in config[entry] and entry not in ['database', 'genologics']:
          unmade_fldr = config[entry]
          if not pathlib.Path(unmade_fldr).exists():
            os.makedirs(unmade_fldr)
            logger.info("Created path {}".format(unmade_fldr))

        #level two
        elif isinstance(config[entry], collections.Mapping):
          for thing in config[entry].keys():
            if isinstance(config[entry][thing], str) and '/' in config[entry][thing] and entry not in ['database', 'genologics']:
              unmade_fldr = config[entry][thing]
              if not pathlib.Path(unmade_fldr).exists():
                os.makedirs(unmade_fldr)
                logger.info("Created path {}".format(unmade_fldr))

  except Exception as e:
    print("Config error: {}".format(str(e)))
    pass
