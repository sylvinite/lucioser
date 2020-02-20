#!/usr/bin/env python

import builtins
import click
import json
import pathlib
import pdb
import pytest
import mock
import os

from microSALT import __version__

from click.testing import CliRunner
from distutils.sysconfig import get_python_lib
from unittest.mock import patch, mock_open

from microSALT import preset_config, logger
from microSALT.cli import root

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
  testdata = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata.json'))
  #Check if release install exists
  for entry in os.listdir(get_python_lib()):
    if 'microSALT-' in entry:
      testdata = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'testdata/testdata.json'))
  return testdata

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
@patch('microSALT.cli.os.path.isdir')
@patch('microSALT.utils.job_creator.glob.glob')
def test_analyse(jc_glob, isdir, gzip, listdir, subproc, runner, config, path_testdata):
  #Sets up subprocess mocking
  process_mock = mock.Mock()
  attrs = {'communicate.return_value': ('output 123456789', 'error')}
  process_mock.configure_mock(**attrs)
  subproc.return_value = process_mock

  jc_glob.return_value = ['AAA1234A1','AAA1234A2']

  isdir.return_value = True
  listdir.return_value = ["ACC6438A3_HVMHWDSXX_L1_1.fastq.gz", "ACC6438A3_HVMHWDSXX_L1_2.fastq.gz", "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz", "ACC6438A3_HVMHWDSXX_L2_2.fastq.gz"]

  #All subcommands
  base_invoke = runner.invoke(root, ['analyse'])
  assert base_invoke.exit_code == 2
  #Exhaustive parameter test
  typical_run = runner.invoke(root, ['analyse', path_testdata, '--input', '/tmp/AAA1234', '--config', config, '--email', '2@2.com'])
  assert typical_run.exit_code == 0
  dry_run = runner.invoke(root, ['analyse', path_testdata, '--input', '/tmp/AAA1234', '--dry'])
  assert dry_run.exit_code == 0
  special_run_settings = {'trimmed':False, 'careful':False,'skip_update':True}
  special_run = runner.invoke(root, ['analyse', path_testdata, '--track', 'qc', '--skip_update', '--untrimmed', '--uncareful', '--input', '/tmp/AAA1234'])
  assert special_run.exit_code == 0


@patch('microSALT.utils.job_creator.Job_Creator.create_project')
@patch('microSALT.utils.reporter.Reporter.start_web')
@patch('multiprocessing.Process.terminate')
@patch('multiprocessing.Process.join')
@patch('microSALT.utils.reporter.requests.get')
@patch('microSALT.utils.reporter.smtplib')
@patch('os.listdir')
@patch('os.path.isdir')
def test_finish(isdir, listdir, smtp, reqs_get, proc_join, proc_term, webstart, create_projct, runner, config, path_testdata):
  listdir.return_value = ['AAA1234A1', 'AAA1234A2' , 'AAA1234A3']
  isdir.return_value = True

  #All subcommands
  base_invoke = runner.invoke(root, ['utils', 'finish'])
  assert base_invoke.exit_code == 2
   #Exhaustive parameter test
  typical_run = runner.invoke(root, ['utils', 'finish', path_testdata, '--email', '2@2.com', '--input', '/tmp/AAA1234_2019.8.12_11.25.2', '--config', config, '--report', 'default'])
  assert typical_run.exit_code == 0
  special_run = runner.invoke(root, ['utils', 'finish', path_testdata, '--report', 'qc'])
  assert special_run.exit_code == 0
  unique_report = runner.invoke(root, ['utils', 'finish', path_testdata, '--report', 'motif_overview'])
  assert unique_report.exit_code == 0


@patch('microSALT.utils.reporter.Reporter.start_web')
@patch('multiprocessing.Process.terminate')
@patch('multiprocessing.Process.join')
@patch('microSALT.utils.reporter.requests.get')
@patch('microSALT.utils.reporter.smtplib')
def test_report(smtplib, reqget, join, term, webstart, runner, path_testdata):
  base_invoke = runner.invoke(root, ['utils', 'report'])
  assert base_invoke.exit_code == 2

  #Exhaustive parameter test
  for rep_type in ['default','typing','motif_overview','qc','json_dump','st_update']:
    normal_report = runner.invoke(root, ['utils', 'report', path_testdata,'--type', rep_type, '--email', '2@2.com', '--output', '/tmp/'])
    assert normal_report.exit_code == 0
    collection_report = runner.invoke(root, ['utils', 'report',path_testdata, '--type', rep_type, '--collection'])
    assert collection_report.exit_code == 0

@patch('microSALT.utils.reporter.Reporter.start_web')
@patch('multiprocessing.Process.terminate')
@patch('multiprocessing.Process.join')
@patch('microSALT.utils.reporter.requests.get')
@patch('microSALT.utils.reporter.smtplib')
def test_resync(smtplib, reqget, join, term, webstart, runner):
  a = runner.invoke(root, ['utils', 'resync', 'overwrite', 'AAA1234A1'])
  assert a.exit_code == 0
  b = runner.invoke(root, ['utils', 'resync', 'overwrite', 'AAA1234A1', '--force'])
  assert b.exit_code == 0

  #Exhaustive parameter test
  for rep_type in ['list', 'report']:
    typical_work = runner.invoke(root, ['utils', 'resync', 'review', '--email', '2@2.com', '--type', rep_type])
    assert typical_work.exit_code == 0
    delimited_work = runner.invoke(root, ['utils', 'resync', 'review', '--skip_update', '--customer', 'custX', '--type', rep_type])
    assert delimited_work.exit_code == 0

def test_refer(runner):
  list_invoke = runner.invoke(root, ['utils', 'refer', 'list'])
  assert list_invoke.exit_code == 0

  a = runner.invoke(root, ['utils', 'refer', 'add', 'Homosapiens_Trams'])
  assert a.exit_code == 0
  b = runner.invoke(root, ['utils', 'refer', 'add', 'Homosapiens_Trams', '--force'])
  assert b.exit_code == 0

@patch('microSALT.utils.reporter.Reporter.start_web')
def test_view(webstart, runner):
  view = runner.invoke(root, ['utils', 'view'])
  assert view.exit_code == 0

@patch('subprocess.Popen')
def test_autobatch(subproc, runner):
  #Sets up subprocess mocking
  process_mock = mock.Mock()
  attrs = {'communicate.return_value': ('"AAA1000_job"', 'error')}
  process_mock.configure_mock(**attrs)
  subproc.return_value = process_mock

  ab = runner.invoke(root, ['utils', 'autobatch', '--dry'])
  assert ab.exit_code == 0

#@patch('os.path.isdir')
#def test_generate(isdir, runner):
#  gent = runner.invoke(root, ['utils', 'generate', '--input', '/tmp/AAA1234'])
#  assert gent.exit_code == 0
