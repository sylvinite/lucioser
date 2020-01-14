#!/usr/bin/env python

import click
import pytest
import mock
import os

from microSALT import __version__

from click.testing import CliRunner
from unittest.mock import patch

from microSALT import config, logger
from microSALT.cli import root

@pytest.fixture
def runner():
  runnah = CliRunner()
  return runnah

@pytest.fixture
def config():
  config = "{}/configExample.json".format(os.path.dirname(os.getcwd()))
  return config

def test_version(runner):
  res = runner.invoke(root, '--version')
  assert res.exit_code == 0
  assert __version__ in res.stdout

@patch('microSALT.store.lims_fetcher.Lims.check_version')
def test_groups(check_version, runner):
  """These groups should only return the help text"""
  base = runner.invoke(root, ['analyse'])
  assert base.exit_code == 0
  base = runner.invoke(root, ['utils'])
  assert base.exit_code == 0
  base_invoke = runner.invoke(root, ['utils', 'resync'])
  assert base_invoke.exit_code == 0
  base_invoke = runner.invoke(root, ['utils', 'refer'])
  assert base_invoke.exit_code == 0

@patch('subprocess.Popen')
@patch('os.listdir')
@patch('os.path.isdir')
@patch('microSALT.store.lims_fetcher.Lims.get_samples')
@patch('microSALT.store.lims_fetcher.Lims.check_version')
def test_analyse(check_version, get_samples, isdir, listdir, subproc, runner, config):
  #Sets up subprocess mocking
  process_mock = mock.Mock()
  attrs = {'communicate.return_value': ('output', 'error')}
  process_mock.configure_mock(**attrs)
  subproc.return_value = process_mock 

  #All subcommands
  for analysis_type in ['sample', 'project', 'collection']:
    base_invoke = runner.invoke(root, ['analyse', analysis_type])
    assert base_invoke.exit_code == 2

    #Exhaustive parameter test
    typical_run = runner.invoke(root, ['analyse', analysis_type, 'AAA1234', '--input', '/tmp/AAA1234', '--config', config, '--email', '2@2.com'])
    assert typical_run.exit_code == 0
    dry_run = runner.invoke(root, ['analyse', analysis_type, 'AAA1234', '--dry'])
    assert dry_run.exit_code == 0
    special_run = runner.invoke(root, ['analyse', analysis_type, 'AAA1234', '--qc_only', '--skip_update', '--untrimmed', '--uncareful'])
    assert special_run.exit_code == 0

@patch('microSALT.utils.reporter.Reporter.start_web')
@patch('microSALT.store.lims_fetcher.Lims.check_version')
def test_finish(check_version, webstart, runner, config):
  #All subcommands
  for analysis_type in ['sample', 'project', 'collection']:
    base_invoke = runner.invoke(root, ['utils', 'finish', analysis_type])
    assert base_invoke.exit_code == 2

    #Exhaustive parameter test
    typical_run = runner.invoke(root, ['utils', 'finish', analysis_type, 'AAA1234', '--email', '2@2.com', '--input', '/tmp/', '--config', config, '--report', 'default'])
    assert typical_run.exit_code == -1
    special_run = runner.invoke(root, ['utils', 'finish', analysis_type, 'AAA1234', '--rerun', '--report', 'qc'])
    assert special_run.exit_code == -1
    if analysis_type == 'collection':
      unique_report = runner.invoke(root, ['utils', 'finish', analysis_type, 'AAA1234', '--report', 'motif_overview'])
      assert unique_report.exit_code == -1


@patch('microSALT.utils.reporter.Reporter.start_web')
@patch('multiprocessing.Process.terminate')
@patch('multiprocessing.Process.join')
@patch('microSALT.utils.reporter.requests.get')
@patch('microSALT.utils.reporter.smtplib')
@patch('microSALT.utils.reporter.LIMS_Fetcher')
@patch('microSALT.store.lims_fetcher.Lims.check_version')
def test_report(check_version, LF, smtplib, reqget, join, term, webstart, runner):
  base_invoke = runner.invoke(root, ['utils', 'report'])
  assert base_invoke.exit_code == 2

  #Exhaustive parameter test
  for rep_type in ['default','typing','motif_overview','qc','json_dump','st_update']:
    normal_report = runner.invoke(root, ['utils', 'report', 'AAA1234', '--type', rep_type, '--email', '2@2.com', '--output', '/tmp/'])
    assert normal_report.exit_code == 0
    collection_report = runner.invoke(root, ['utils', 'report', 'AAA1234', '--type', rep_type, '--collection'])
    assert collection_report.exit_code == 0

@patch('microSALT.utils.reporter.Reporter.start_web')
@patch('multiprocessing.Process.terminate')
@patch('multiprocessing.Process.join')
@patch('microSALT.utils.reporter.requests.get')
@patch('microSALT.utils.reporter.smtplib')
@patch('microSALT.store.lims_fetcher.Lims.check_version')
def test_resync(check_version, smtplib, reqget, join, term, webstart, runner):
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

@patch('microSALT.store.lims_fetcher.Lims.check_version')
def test_refer(check_version, runner):
  list_invoke = runner.invoke(root, ['utils', 'refer', 'list'])
  assert list_invoke.exit_code == 0

  a = runner.invoke(root, ['utils', 'refer', 'add', 'Homosapiens_Trams'])
  assert a.exit_code == 0
  b = runner.invoke(root, ['utils', 'refer', 'add', 'Homosapiens_Trams', '--force'])
  assert b.exit_code == 0

@patch('microSALT.utils.reporter.Reporter.start_web')
@patch('microSALT.store.lims_fetcher.Lims.check_version')
def test_view(check_version, webstart, runner):
  view = runner.invoke(root, ['utils', 'view'])
  assert view.exit_code == 0
