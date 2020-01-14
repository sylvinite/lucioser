#!/usr/bin/env python

import collections
import os
import pathlib
import pytest

from microSALT import config

@pytest.fixture
def exp_config():
  precon = \
  {
    'slurm_header': 
      {'time','threads', 'qos', 'job_prefix','project', 'type'},
    'regex':
      {'file_pattern', 'mail_recipient', 'verified_organisms'},
    'folders':
      {'results', 'reports', 'log_file', 'seqdata', 'profiles', 'references', 'resistances', 'genomes', 'expec', 'adapters', 'cgmlst'},
    'threshold':
      {'mlst_id', 'mlst_novel_id', 'mlst_span', 'cgmlst_id', 'cgmlst_coverage', 'motif_id', 'motif_span', 'total_reads_warn', 'total_reads_fail', 'NTC_total_reads_warn', \
                       'NTC_total_reads_fail', 'mapped_rate_warn', 'mapped_rate_fail', 'duplication_rate_warn', 'duplication_rate_fail', 'insert_size_warn', 'insert_size_fail', \
                       'average_coverage_warn', 'average_coverage_fail', 'bp_10x_warn', 'bp_10x_fail', 'bp_30x_warn', 'bp_50x_warn', 'bp_100x_warn'},
    'database':
      {'SQLALCHEMY_DATABASE_URI' ,'SQLALCHEMY_TRACK_MODIFICATIONS' , 'DEBUG'},
    'genologics':
      {'baseuri', 'username', 'password'},
    'dry': True,
  }
  return precon

def test_existence(exp_config):
  """Checks that the configuration contains certain key variables"""

  #level one
  config_level_one = config.keys()
  for entry in exp_config.keys():
    if entry != 'dry':
      assert entry in config_level_one

      #level two
      if isinstance(config[entry], collections.Mapping):
        config_level_two = config[entry].keys()
        for thing in exp_config[entry]:
          assert thing in config_level_two

def test_reverse_existence(exp_config):
  """Check that the configuration doesnt contain outdated variables"""

  #level one
  config_level_one = exp_config.keys()
  for entry in config.keys():
    if entry not in ['_comment', 'rerun']:
      assert entry in config_level_one

      #level two
      config_level_two = exp_config[entry]
      if isinstance(config[entry], collections.Mapping):
        for thing in config[entry].keys():
          if thing != '_comment':
            assert thing in config_level_two

#def test_type(exp_config):
#  """Verify that each variable uses the correct format"""
#  pass

def test_paths(exp_config):
  """Tests existence for all paths mentioned in variables"""
  #level one
  for entry in config.keys():
    if entry != '_comment':
      if isinstance(config[entry], str) and '/' in config[entry] and entry not in ['database', 'genologics']:
        unmade_fldr = config[entry][thing]
        assert (pathlib.Path(unmade_fldr).exists())
    
      #level two
      elif isinstance(config[entry], collections.Mapping):
        for thing in config[entry].keys():
          if isinstance(config[entry][thing], str) and '/' in config[entry][thing] and entry not in ['database', 'genologics']:
            unmade_fldr = config[entry][thing]
            assert (pathlib.Path(unmade_fldr).exists())
    
