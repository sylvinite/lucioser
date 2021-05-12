#!/usr/bin/env python

import json
import os
import pathlib
import pdb
import pytest
import requests
import re
import runpy
import time

from distutils.sysconfig import get_python_lib
from unittest.mock import patch

from microSALT.utils.reporter import Reporter
from microSALT import preset_config, logger
from microSALT.cli import root
from microSALT.server.views import *
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
def testdata():
  testdata = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata/sampleinfo_samples.json'))
  #Check if release install exists
  for entry in os.listdir(get_python_lib()):
    if 'microSALT-' in entry:
      testdata = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'testdata/sampleinfo_samples.json'))
  with open(testdata) as json_file:
    data = json.load(json_file)
  return data

@pytest.fixture
def report_obj(testdata):
  report = Reporter(config=preset_config, log=logger, sampleinfo=testdata)
  return report

@pytest.fixture
def appscript():
  script = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'microSALT/server/app.py'))
  return script

def test_webserver(report_obj):
  report_obj.start_web()
  report_obj.kill_flask()

def test_appobject(appscript):
  runpy.run_path(appscript)   

def test_pages(report_obj, mock_db):
  report_obj.start_web()
  #Valid pages with available data
  time.sleep(3)
  a = requests.get("http://127.0.0.1:5000/", allow_redirects=True)
  assert a.status_code == 200
  
  time.sleep(0.15)
  b = requests.get("http://127.0.0.1:5000/microSALT/", allow_redirects=True)
  assert b.status_code == 200

  time.sleep(0.15)
  c = requests.get("http://127.0.0.1:5000/microSALT/AAA1234", allow_redirects=True)
  assert c.status_code == 200

  time.sleep(0.15)
  e = requests.get("http://127.0.0.1:5000/microSALT/AAA1234/typing/all", allow_redirects=True)
  assert e.status_code in [200, 500]

  time.sleep(0.15)
  d = requests.get("http://127.0.0.1:5000/microSALT/AAA1234/qc", allow_redirects=True)
  assert d.status_code in [200, 500]

  #Valid pages with unavailable data
  time.sleep(0.15)
  f = requests.get("http://127.0.0.1:5000/microSALT/AAA1234/typing/escherichia_coli", allow_redirects=True)
  assert f.status_code in [200, 500]

  time.sleep(0.15)
  g = requests.get("http://127.0.0.1:5000/microSALT/STtracker/all", allow_redirects=True)
  assert g.status_code in [200, 500]

  time.sleep(0.15)
  h = requests.get("http://127.0.0.1:5000/microSALT/STtracker/cust000", allow_redirects=True)
  assert h.status_code in [200, 500]

  report_obj.kill_flask()

@patch('microSALT.server.views.render_template')
def test_index_views(renderpatch, mock_db):
  renderpatch.return_value = "ok"
  start = start_page()
  assert start == "ok"
  reroute = reroute_page()
  assert reroute == "ok"

@patch('microSALT.server.views.render_template')
def test_project_views(renderpatch, mock_db):
  renderpatch.return_value = "ok"
  a = project_page("AAA1234")
  assert a == "ok"
  b = alignment_page("AAA1234")
  assert b == "ok"
  c = typing_page("AAA1234","all")
  assert c == "ok"

@patch('microSALT.server.views.gen_add_info')
@patch('microSALT.server.views.render_template')
def test_tracker_view(renderpatch, addinfo, mock_db):
  renderpatch.return_value = "ok"
  a = STtracker_page("cust000")
  assert a == "ok"
   
