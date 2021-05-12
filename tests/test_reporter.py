#!/usr/bin/env python

import datetime
import glob
import json
import logging
import os
import pathlib
import pdb
import pytest
import re
import sys

from distutils.sysconfig import get_python_lib
from unittest.mock import patch

from microSALT import preset_config, logger
from microSALT.utils.reporter import Reporter
from microSALT.utils.referencer import Referencer
from microSALT.store.db_manipulator import DB_Manipulator

def unpack_db_json(filename):
  testdata = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata/{}'.format(filename)))
  #Check if release install exists
  for entry in os.listdir(get_python_lib()):
    if 'microSALT-' in entry:
      testdata = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'testdata/{}'.format(filename)))
  with open(testdata) as json_file:
    data = json.load(json_file)
  return data

@pytest.fixture
def mock_db():
  db_file = re.search('sqlite:///(.+)', preset_config['database']['SQLALCHEMY_DATABASE_URI']).group(1)
  dbm = DB_Manipulator(config=preset_config,log=logger)
  dbm.create_tables()

  for antry in unpack_db_json('sampleinfo_projects.json'):
    dbm.add_rec(antry, 'Projects')
  for entry in unpack_db_json('sampleinfo_mlst.json'):
    dbm.add_rec(entry, 'Seq_types')
  for bentry in unpack_db_json('sampleinfo_resistance.json'):
    dbm.add_rec(bentry, 'Resistances')
  for centry in unpack_db_json('sampleinfo_custom.json'):
    dbm.add_rec(centry, 'Custom_targets')
  for dentry in unpack_db_json('sampleinfo_reports.json'):
    dbm.add_rec(dentry, 'Reports')
  return dbm

@pytest.fixture
def reporter():
  reporter_obj = Reporter(config=preset_config, log=logger, sampleinfo=unpack_db_json('sampleinfo_samples.json')[0], name="MIC1234A1", output="/tmp/MLST")
  return reporter_obj

def test_motif(mock_db, reporter):
  reporter.create_subfolders()
  reporter.gen_motif(motif="resistance")
  assert len( glob.glob("{}/AAA1234_resistance*".format(reporter.output))) > 0

  reporter.create_subfolders()
  reporter.gen_motif(motif="custom")
  assert len( glob.glob("{}/AAA1234_custom*".format(reporter.output))) > 0


def test_deliveryreport(mock_db, reporter):
  reporter.create_subfolders()
  reporter.gen_delivery()
  assert len( glob.glob("{}/deliverables/999999_deliverables.yaml".format(preset_config['folders']['reports']))) > 0

def test_jsonreport(mock_db, reporter):
  reporter.create_subfolders()
  reporter.gen_json()
  assert len( glob.glob("{}/json/AAA1234.json".format(preset_config['folders']['reports']))) > 0

def test_gen_qc(mock_db, reporter):
  reporter.name = "name_that_do_not_exist"
  with pytest.raises(Exception):
    reporter.gen_qc()

def test_gen_typing(mock_db, reporter):
  reporter.name = "name_that_do_not_exist"
  with pytest.raises(Exception):
    reporter.gen_typing()

def test_gen_motif(caplog, reporter):
  caplog.clear()
  reporter.gen_motif(motif="unrecognized")
  assert "Invalid motif type" in caplog.text
  caplog.clear()
  reporter.output = "/path/that/do/not/exists/"
  reporter.gen_motif()
  assert "Gen_motif unable to produce" in caplog.text

def test_gen_json(caplog, reporter):
  caplog.clear()
  reporter.output = "/path/that/do/not/exists/"
  preset_config["folders"]["reports"] = "/path/that/do/not/exists/"
  reporter.config = preset_config
  reporter.gen_json()
  assert "Gen_json unable to produce" in caplog.text

def test_report(caplog, reporter):
  caplog.clear()
  reporter.type = "type_not_mentioned_in_list"
  with pytest.raises(Exception):
    reporter.report()
    assert "Report function recieved invalid format" in caplog.text

@patch('microSALT.utils.reporter.Reporter.start_web')
def test_restart_web(sw, reporter):
  reporter.restart_web()

def test_constructor():
  sample_info = [
                {"CG_ID_project" : "AAA1234","CG_ID_sample" : "AAA1234A1","Customer_ID_project" :
                  "999999","Customer_ID_sample" : "XXX0000Y1","Customer_ID" : "cust000","application_tag" :
                  "NONE","date_arrival" : "0001-01-01 00:00:00","date_libprep" : "0001-01-01 00:00:00",
                  "date_sequencing" : "0001-01-01 00:00:00","method_libprep" : "Not in LIMS",
                  "method_sequencing" : "Not in LIMS","organism" : "Staphylococcus aureus","priority" :
                  "standard","reference" : "AP017922.1"},
                {"CG_ID_project" : "Something_else_than_AAA1234","CG_ID_sample" : "AAA1234A2",
                 "Customer_ID_project" : "999999","Customer_ID_sample" : "XXX0000Y1","Customer_ID" :
                 "cust000","application_tag" : "NONE","date_arrival" : "0001-01-01 00:00:00","date_libprep" :
                 "0001-01-01 00:00:00","date_sequencing" : "0001-01-01 00:00:00","method_libprep" :
                 "Not in LIMS","method_sequencing" : "Not in LIMS","organism" : "Escherichia coli",
                 "priority" : "standard","reference" : "NC_011751.1"}
                ]
  reporter_obj = Reporter(config=preset_config, log=logger, sampleinfo=sample_info,
                          name="MIC1234A1", output="/tmp/MLST")
