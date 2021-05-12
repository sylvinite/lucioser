#!/usr/bin/env python

import builtins
import click
import json
import logging
import pathlib
import pdb
import pytest
import re
import mock
import os
import sys

from microSALT import __version__

from click.testing import CliRunner
from distutils.sysconfig import get_python_lib
from unittest.mock import patch, mock_open

from microSALT import preset_config, logger
from microSALT.cli import root
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
def dbm():
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


@pytest.fixture(autouse=True)
def no_requests(monkeypatch):
    """Remove requests.sessions.Session.request for all tests."""
    monkeypatch.delattr("requests.sessions.Session.request")

@pytest.fixture
def runner():
  runnah = CliRunner()
  return runnah

@pytest.fixture
def config():
  config = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'configExample.json'))
  #Check if release install exists
  for entry in os.listdir(get_python_lib()):
    if 'microSALT-' in entry:
      config = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'testdata/configExample.json'))
  return config

@pytest.fixture
def path_testdata():
  testdata = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata/sampleinfo_samples.json'))
  #Check if release install exists
  for entry in os.listdir(get_python_lib()):
    if 'microSALT-' in entry:
      testdata = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'testdata/sampleinfo_samples.json'))
  return testdata

@pytest.fixture
def path_testproject():
  testproject = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata/AAA1234_2000.1.2_3.4.5'))
  #Check if release install exists
  for entry in os.listdir(get_python_lib()):
    if 'microSALT-' in entry:
      testproject = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'testproject/AAA1234_2000.1.2_3.4.5'))
  return testproject


def test_version(runner):
  res = runner.invoke(root, '--version')
  assert res.exit_code == 0
  assert __version__ in res.stdout

def test_groups(runner):
  """These groups should only return the help text"""
  base = runner.invoke(root, ['utils'])
  assert base.exit_code == 0
  base_invoke = runner.invoke(root, ['utils', 'resync'])
  assert base_invoke.exit_code == 0
  base_invoke = runner.invoke(root, ['utils', 'refer'])
  assert base_invoke.exit_code == 0


@patch('subprocess.Popen')
@patch('os.listdir')
@patch('gzip.open')
@patch('microSALT.utils.job_creator.glob.glob')
@patch('microSALT.cli.os.path.isdir')
def test_analyse_dry(isdir, jc_glob, gzip, listdir, subproc, runner, config, path_testdata, caplog, dbm):
  caplog.clear()
  caplog.set_level(logging.DEBUG, logger="main_logger")

  #Sets up subprocess mocking
  process_mock = mock.Mock()
  attrs = {'communicate.return_value': ('output 123456789', 'error')}
  process_mock.configure_mock(**attrs)
  subproc.return_value = process_mock
  isdir.return_value = True

  jc_glob.return_value = ['AAA1234A1','AAA1234A2']

#  isdir.return_value = True
  listdir.return_value = ["ACC6438A3_HVMHWDSXX_L1_1.fastq.gz", "ACC6438A3_HVMHWDSXX_L1_2.fastq.gz", "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz", "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz"]

  #All subcommands
  base_invoke = runner.invoke(root, ['analyse'])
  assert base_invoke.exit_code == 2
  #Exhaustive parameter test
  dry_run = runner.invoke(root, ['analyse', path_testdata, '--input', '/tmp/AAA1234', '--dry'])
  assert dry_run.exit_code == 0
  assert "INFO - Execution finished!" in caplog.text
  caplog.clear()

@patch('subprocess.Popen')
@patch('os.listdir')
@patch('gzip.open')
@patch('microSALT.utils.job_creator.glob.glob')
@patch('microSALT.cli.os.path.isdir')
def test_analyse_typical(isdir, jc_glob, gzip, listdir, subproc, runner, config, path_testdata, caplog, dbm):
  caplog.clear()
  caplog.set_level(logging.DEBUG, logger="main_logger")

  #Sets up subprocess mocking
  process_mock = mock.Mock()
  attrs = {'communicate.return_value': ('output 123456789', 'error')}
  process_mock.configure_mock(**attrs)
  subproc.return_value = process_mock
  isdir.return_value = True

  jc_glob.return_value = ['AAA1234A1','AAA1234A2']

#  isdir.return_value = True
  listdir.return_value = ["ACC6438A3_HVMHWDSXX_L1_1.fastq.gz", "ACC6438A3_HVMHWDSXX_L1_2.fastq.gz", "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz", "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz"]

  typical_run = runner.invoke(root, ['analyse', path_testdata, '--input', '/tmp/AAA1234', '--config', config, '--email', '2@2.com'])
  assert typical_run.exit_code == 0
  assert "INFO - Execution finished!" in caplog.text
  caplog.clear()

@patch('subprocess.Popen')
@patch('os.listdir')
@patch('gzip.open')
@patch('microSALT.utils.job_creator.glob.glob')
@patch('microSALT.cli.os.path.isdir')
def test_analyse_exhaustive(isdir, jc_glob, gzip, listdir, subproc, runner, config, path_testdata, caplog, dbm):
  caplog.clear()
  caplog.set_level(logging.DEBUG, logger="main_logger")

  #Sets up subprocess mocking
  process_mock = mock.Mock()
  attrs = {'communicate.return_value': ('output 123456789', 'error')}
  process_mock.configure_mock(**attrs)
  subproc.return_value = process_mock
  isdir.return_value = True

  jc_glob.return_value = ['AAA1234A1','AAA1234A2']

#  isdir.return_value = True
  listdir.return_value = ["ACC6438A3_HVMHWDSXX_L1_1.fastq.gz", "ACC6438A3_HVMHWDSXX_L1_2.fastq.gz", "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz", "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz"]

  special_run_settings = {'trimmed':False, 'careful':False,'skip_update':True}
  special_run = runner.invoke(root, ['analyse', path_testdata, '--skip_update', '--untrimmed', '--uncareful', '--input', '/tmp/AAA1234'])
  assert special_run.exit_code == 0
  assert "INFO - Execution finished!" in caplog.text

@patch('microSALT.utils.job_creator.Job_Creator.create_project')
@patch('microSALT.utils.reporter.Reporter.start_web')
@patch('multiprocessing.Process.terminate')
@patch('multiprocessing.Process.join')
@patch('microSALT.utils.reporter.requests.get')
@patch('microSALT.utils.reporter.smtplib')
@patch('microSALT.cli.os.path.isdir')
def test_finish_typical(isdir, smtp, reqs_get, proc_join, proc_term, webstart, create_projct, runner, config, path_testdata, path_testproject, caplog, dbm):
  caplog.set_level(logging.DEBUG, logger="main_logger")
  caplog.clear()

  isdir.return_value = True

  #All subcommands
  base_invoke = runner.invoke(root, ['utils', 'finish'])
  assert base_invoke.exit_code == 2
   #Exhaustive parameter test
  typical_run = runner.invoke(root, ['utils', 'finish', path_testdata, '--email', '2@2.com', '--config', config, '--report', 'default', '--output', '/tmp/', '--input', path_testproject])
  assert typical_run.exit_code == 0
  assert "INFO - Execution finished!" in caplog.text
  caplog.clear()

@patch('microSALT.utils.job_creator.Job_Creator.create_project')
@patch('microSALT.utils.reporter.Reporter.start_web')
@patch('multiprocessing.Process.terminate')
@patch('multiprocessing.Process.join')
@patch('microSALT.utils.reporter.requests.get')
@patch('microSALT.utils.reporter.smtplib')
@patch('microSALT.cli.os.path.isdir')
def test_finish_qc(isdir, smtp, reqs_get, proc_join, proc_term, webstart, create_projct, runner, config, path_testdata, path_testproject, caplog, dbm):
  caplog.set_level(logging.DEBUG, logger="main_logger")
  caplog.clear()

  isdir.return_value = True

  special_run = runner.invoke(root, ['utils', 'finish', path_testdata, '--report', 'qc', '--output', '/tmp/', '--input', path_testproject])
  assert special_run.exit_code == 0
  assert "INFO - Execution finished!" in caplog.text
  caplog.clear()

@patch('microSALT.utils.job_creator.Job_Creator.create_project')
@patch('microSALT.utils.reporter.Reporter.start_web')
@patch('multiprocessing.Process.terminate')
@patch('multiprocessing.Process.join')
@patch('microSALT.utils.reporter.requests.get')
@patch('microSALT.utils.reporter.smtplib')
@patch('microSALT.cli.os.path.isdir')
def test_finish_motif(isdir, smtp, reqs_get, proc_join, proc_term, webstart, create_projct, runner, config, path_testdata, path_testproject, caplog, dbm):
  caplog.set_level(logging.DEBUG, logger="main_logger")
  caplog.clear()

  isdir.return_value = True

  unique_report = runner.invoke(root, ['utils', 'finish', path_testdata, '--report', 'motif_overview', '--output', '/tmp/', '--input', path_testproject])
  assert unique_report.exit_code == 0
  assert "INFO - Execution finished!" in caplog.text
  caplog.clear()


@patch('microSALT.utils.reporter.Reporter.start_web')
@patch('multiprocessing.Process.terminate')
@patch('multiprocessing.Process.join')
@patch('microSALT.utils.reporter.requests.get')
@patch('microSALT.utils.reporter.smtplib')
def test_report(smtplib, reqget, join, term, webstart, runner, path_testdata, caplog, dbm):
  caplog.set_level(logging.DEBUG, logger="main_logger")
  caplog.clear()

  base_invoke = runner.invoke(root, ['utils', 'report'])
  assert base_invoke.exit_code == 2

  #Exhaustive parameter test
  for rep_type in ['default','typing','motif_overview','qc','json_dump','st_update']:
    normal_report = runner.invoke(root, ['utils', 'report', path_testdata,'--type', rep_type, '--email', '2@2.com', '--output', '/tmp/'])
    assert normal_report.exit_code == 0
    assert "INFO - Execution finished!" in caplog.text
    caplog.clear()
    collection_report = runner.invoke(root, ['utils', 'report',path_testdata, '--type', rep_type, '--collection', '--output', '/tmp/'])
    assert collection_report.exit_code == 0
    assert "INFO - Execution finished!" in caplog.text
    caplog.clear()


@patch('microSALT.utils.reporter.Reporter.start_web')
@patch('multiprocessing.Process.terminate')
@patch('multiprocessing.Process.join')
@patch('microSALT.utils.reporter.requests.get')
@patch('microSALT.utils.reporter.smtplib')
def test_resync_overwrite(smtplib, reqget, join, term, webstart, runner, caplog, dbm):
  caplog.set_level(logging.DEBUG, logger="main_logger")
  caplog.clear()

  a = runner.invoke(root, ['utils', 'resync', 'overwrite', 'AAA1234A1'])
  assert a.exit_code == 0
  assert "INFO - Execution finished!" in caplog.text
  caplog.clear()
  b = runner.invoke(root, ['utils', 'resync', 'overwrite', 'AAA1234A1', '--force'])
  assert b.exit_code == 0
  assert "INFO - Execution finished!" in caplog.text
  caplog.clear()

@patch('microSALT.utils.reporter.Reporter.start_web')
@patch('multiprocessing.Process.terminate')
@patch('multiprocessing.Process.join')
@patch('microSALT.utils.reporter.requests.get')
@patch('microSALT.utils.reporter.smtplib')
def test_resync_review(smtplib, reqget, join, term, webstart, runner, caplog, dbm):
  caplog.set_level(logging.DEBUG, logger="main_logger")
  caplog.clear()

  #Exhaustive parameter test
  for rep_type in ['list', 'report']:
    typical_work = runner.invoke(root, ['utils', 'resync', 'review' ,'--email', '2@2.com', '--type', rep_type, '--output', '/tmp/'])
    assert typical_work.exit_code == 0
    assert "INFO - Execution finished!" in caplog.text
    caplog.clear()
    delimited_work = runner.invoke(root, ['utils', 'resync', 'review', '--skip_update', '--customer', 'custX', '--type', rep_type, '--output', '/tmp/'])
    assert delimited_work.exit_code == 0
    assert "INFO - Execution finished!" in caplog.text
    caplog.clear()

def test_refer(runner, caplog, dbm):
  caplog.set_level(logging.DEBUG, logger="main_logger")

  list_invoke = runner.invoke(root, ['utils', 'refer', 'observe'])
  assert list_invoke.exit_code == 0

  a = runner.invoke(root, ['utils', 'refer', 'add', 'Homosapiens_Trams'])
  assert a.exit_code == 0
  #assert "INFO - Execution finished!" in caplog.text
  caplog.clear()
  b = runner.invoke(root, ['utils', 'refer', 'add', 'Homosapiens_Trams', '--force'])
  assert b.exit_code == 0
  #assert "INFO - Execution finished!" in caplog.text
  caplog.clear()

@patch('microSALT.utils.reporter.Reporter.start_web')
def test_view(webstart, runner, caplog, dbm):
  caplog.set_level(logging.DEBUG, logger="main_logger")

  view = runner.invoke(root, ['utils', 'view'])
  assert view.exit_code == 0
  #assert "INFO - Execution finished!" in caplog.text
  caplog.clear()

@patch('os.path.isdir')
def test_generate(isdir, runner, caplog, dbm):
  caplog.set_level(logging.DEBUG, logger="main_logger")
  gent = runner.invoke(root, ['utils', 'generate', '--input', '/tmp/'])
  assert gent.exit_code == 0
  fent = runner.invoke(root, ['utils', 'generate'])
  assert fent.exit_code == 0

