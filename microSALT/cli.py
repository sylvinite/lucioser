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
from microSALT import __version__, preset_config, logger, wd
from microSALT.utils.scraper import Scraper
from microSALT.utils.job_creator import Job_Creator
from microSALT.utils.reporter import Reporter
from microSALT.utils.referencer import Referencer

default_sampleinfo = {
    "CG_ID_project": "XXX0000",
    "CG_ID_sample": "XXX0000A1",
    "Customer_ID_project": "100100",
    "Customer_ID_sample": "10XY123456",
    "Customer_ID": "cust000",
    "application_tag": "SOMTIN100",
    "date_arrival": "0001-01-01 00:00:00",
    "date_libprep": "0001-01-01 00:00:00",
    "date_sequencing": "0001-01-01 00:00:00",
    "method_libprep": "Not in LIMS",
    "method_sequencing": "Not in LIMS",
    "organism": "Staphylococcus aureus",
    "priority": "standard",
    "reference": "None",
}


if preset_config == "":
    click.echo(
        "ERROR - No properly set-up config under neither envvar MICROSALT_CONFIG nor ~/.microSALT/config.json. Exiting."
    )
    sys.exit(-1)


def set_cli_config(config):
    if config != "":
        if os.path.exists(config):
            try:
                t = ctx.obj["config"]
                with open(os.path.abspath(config), "r") as conf:
                    ctx.obj["config"] = json.load(conf)
                ctx.obj["config"]["folders"]["adapters"] = t["folders"]["adapters"]
                ctx.obj["config"]["config_path"] = os.path.abspath(config)
            except Exception as e:
                pass


def done():
    click.echo("INFO - Execution finished!")
    logger.debug("INFO - Execution finished!")


def review_sampleinfo(pfile):
    """Reviews sample info. Returns loaded json object"""

    try:
        with open(pfile) as json_file:
            data = json.load(json_file)
    except Exception as e:
        click.echo("Unable to read provided sample info file as json. Exiting..")
        sys.exit(-1)

    if isinstance(data, list):
        for entry in data:
            for k, v in default_sampleinfo.items():
                if not k in entry:
                    click.echo(
                        "WARNING - Parameter {} needs to be provided in sample json. Formatting example: ({})".format(
                            k, v
                        )
                    )
    else:
        for k, v in default_sampleinfo.items():
            if not k in data:
                click.echo(
                    "WARNING - Parameter {} needs to be provided in sample json. Formatting example: ({})".format(
                        k, v
                    )
                )
    return data


@click.group()
@click.version_option(__version__)
@click.pass_context
def root(ctx):
    """microbial Sequence Analysis and Loci-based Typing (microSALT) pipeline """
    ctx.obj = {}
    ctx.obj["config"] = preset_config
    ctx.obj["log"] = logger


@root.command()
@click.argument("sampleinfo_file")
@click.option("--input", help="Full path to input folder", default="")
@click.option("--config", help="microSALT config to override default", default="")
@click.option(
    "--dry",
    help="Builds instance without posting to SLURM",
    default=False,
    is_flag=True,
)
@click.option(
    "--email",
    default=preset_config["regex"]["mail_recipient"],
    help="Forced e-mail recipient",
)
@click.option(
    "--skip_update", default=False, help="Skips downloading of references", is_flag=True
)
@click.option(
    "--force_update", default=False, help="Forces downloading of pubMLST references", is_flag=True
)
@click.option(
    "--untrimmed", help="Use untrimmed input data", default=False, is_flag=True
)
@click.option(
    "--uncareful",
    help="Avoids running SPAdes in careful mode. Sometimes fix assemblies",
    default=False,
    is_flag=True,
)
@click.option("--custom_target", help="Fasta of custom targets to gene check against", default="")
@click.pass_context
def analyse(
    ctx, sampleinfo_file, input, config, dry, email, skip_update, force_update, untrimmed, uncareful, custom_target):
    """Sequence analysis, typing and resistance identification"""
    # Run section
    pool = []
    trimmed = not untrimmed
    careful = not uncareful
    set_cli_config(config)
    ctx.obj["config"]["regex"]["mail_recipient"] = email
    ctx.obj["config"]["dry"] = dry
    if custom_target != "":
        custom_target = os.path.abspath(custom_target)
    if not os.path.isdir(input):
        click.echo("ERROR - Sequence data folder {} does not exist.".format(input))
        ctx.abort()
    for subfolder in os.listdir(input):
        if os.path.isdir("{}/{}".format(input, subfolder)):
            pool.append(subfolder)

    run_settings = {
        "input": input,
        "dry": dry,
        "email": email,
        "skip_update": skip_update,
        "trimmed": not untrimmed,
        "careful": not uncareful,
        "pool": pool,
        "custom_target":custom_target,
    }

    # Samples section
    sampleinfo = review_sampleinfo(sampleinfo_file)
    run_creator = Job_Creator(
        config=ctx.obj["config"],
        log=ctx.obj["log"],
        sampleinfo=sampleinfo,
        run_settings=run_settings
    )

    ext_refs = Referencer(
        config=ctx.obj["config"], log=ctx.obj["log"], sampleinfo=sampleinfo, force=force_update
    )
    click.echo("INFO - Checking versions of references..")
    try:
        if not skip_update:
            ext_refs.identify_new(project=True)
            ext_refs.update_refs()
            click.echo("INFO - Version check done. Creating sbatch jobs")
        else:
            click.echo("INFO - Skipping version check.")
    except Exception as e:
        click.echo("{}".format(e))
    if len(sampleinfo) > 1:
        run_creator.project_job()
    elif len(sampleinfo) == 1:
        run_creator.project_job(single_sample=True)
    else:
        ctx.abort()

    done()


@root.group()
@click.pass_context
def utils(ctx):
    """ Utilities for specific purposes """
    pass


@utils.group()
@click.pass_context
def refer(ctx):
    """ Manipulates MLST organisms """
    pass


@utils.command()
@click.argument("sampleinfo_file")
@click.option("--input", help="Full path to project folder", default="")
@click.option(
    "--track",
    help="Run a specific analysis track",
    default="default",
    type=click.Choice(["default", "typing", "qc", "cgmlst"]),
)
@click.option("--config", help="microSALT config to override default", default="")
@click.option(
    "--dry",
    help="Builds instance without posting to SLURM",
    default=False,
    is_flag=True,
)
@click.option(
    "--email",
    default=preset_config["regex"]["mail_recipient"],
    help="Forced e-mail recipient",
)
@click.option(
    "--skip_update", default=False, help="Skips downloading of references", is_flag=True
)
@click.option(
    "--report",
    default="default",
    type=click.Choice(
        ["default", "typing", "motif_overview", "qc", "json_dump", "st_update", "multiqc"]
    ),
)
@click.option("--output", help="Report output folder", default="")
@click.option("--custom_target", help="Fasta of custom targets to gene check against", default="")
@click.pass_context
def finish(
    ctx, sampleinfo_file, input, track, config, dry, email, skip_update, report, output, custom_target
):
    """Sequence analysis, typing and resistance identification"""
    # Run section
    pool = []
    set_cli_config(config)
    ctx.obj["config"]["regex"]["mail_recipient"] = email
    ctx.obj["config"]["dry"] = dry
    if custom_target != "":
        custom_target =os.path.abspath(custom_target)

    if not os.path.isdir(input):
        click.echo("ERROR - Sequence data folder {} does not exist.".format(input))
        ctx.abort()
    if output == "":
        output = input
    for subfolder in os.listdir(input):
        if os.path.isdir("{}/{}".format(input, subfolder)):
            pool.append(subfolder)

    run_settings = {
        "input": input,
        "track": track,
        "dry": dry,
        "email": email,
        "skip_update": skip_update,
        "custom_target":custom_target,
    }

    # Samples section
    sampleinfo = review_sampleinfo(sampleinfo_file)
    ext_refs = Referencer(
        config=ctx.obj["config"], log=ctx.obj["log"], sampleinfo=sampleinfo
    )
    click.echo("INFO - Checking versions of references..")
    try:
        if not skip_update:
            ext_refs.identify_new(project=True)
            ext_refs.update_refs()
            click.echo("INFO - Version check done. Creating sbatch jobs")
        else:
            click.echo("INFO - Skipping version check.")
    except Exception as e:
        click.echo("{}".format(e))

    res_scraper = Scraper(
        config=ctx.obj["config"], log=ctx.obj["log"], sampleinfo=sampleinfo, run_settings=run_settings, input=input
    )
    if isinstance(sampleinfo, list) and len(sampleinfo) > 1:
        res_scraper.scrape_project()
        # for subfolder in pool:
        #  res_scraper.scrape_sample()
    else:
        res_scraper.scrape_sample()

    codemonkey = Reporter(
        config=ctx.obj["config"],
        log=ctx.obj["log"],
        sampleinfo=sampleinfo,
        output=output,
        collection=True,
    )
    codemonkey.report(report)
    done()


@refer.command()
@click.argument("organism")
@click.option(
    "--force", help="Redownloads existing organism", default=False, is_flag=True
)
@click.pass_context
def add(ctx, organism, force):
    """ Adds a new internal organism from pubMLST """
    referee = Referencer(config=ctx.obj["config"], log=ctx.obj["log"], force=force)
    try:
        referee.add_pubmlst(organism)
    except Exception as e:
        click.echo(e.args[0])
        ctx.abort()
    click.echo("INFO - Checking versions of all references..")
    referee = Referencer(config=ctx.obj["config"], log=ctx.obj["log"], force=force)
    referee.update_refs()


@refer.command()
@click.pass_context
def observe(ctx):
    """ Lists all stored organisms """
    refe = Referencer(config=ctx.obj["config"], log=ctx.obj["log"])
    click.echo("INFO - Currently stored organisms:")
    for org in sorted(refe.existing_organisms()):
        click.echo(org.replace("_", " ").capitalize())


@refer.command()
@click.pass_context
def update(ctx, sampleinfo_file):
    """Updates all ST and Resistance databases"""
    ext_refs = Referencer(config=ctx.obj["config"], log=ctx.obj["log"])
    click.echo("INFO - Checking versions of references..")
    try:
        ext_refs.update_refs()
        click.echo("INFO - Version check done.")
    except Exception as e:
        click.echo("{}".format(e))
    done()

@utils.command()
@click.argument("sampleinfo_file")
@click.option(
    "--email",
    default=preset_config["regex"]["mail_recipient"],
    help="Forced e-mail recipient",
)
@click.option(
    "--type",
    default="default",
    type=click.Choice(
        ["default", "typing", "motif_overview", "qc", "json_dump", "st_update","multiqc"]
    ),
)
@click.option("--output", help="Full path to output folder", default="")
@click.option("--collection", default=False, is_flag=True)
@click.pass_context
def report(ctx, sampleinfo_file, email, type, output, collection):
    """Re-generates report for a project"""
    ctx.obj["config"]["regex"]["mail_recipient"] = email
    sampleinfo = review_sampleinfo(sampleinfo_file)
    codemonkey = Reporter(
        config=ctx.obj["config"],
        log=ctx.obj["log"],
        sampleinfo=sampleinfo,
        output=output,
        collection=collection,
    )
    codemonkey.report(type)
    done()


@utils.command()
@click.pass_context
def view(ctx):
    """Starts an interactive webserver for viewing"""
    codemonkey = Reporter(config=ctx.obj["config"], log=ctx.obj["log"])
    codemonkey.start_web()

@utils.command()
@click.option("--input", help="Full path to project folder", default=os.getcwd())
@click.pass_context
def generate(ctx, input):
    """Creates a blank sample info json for the given input folder"""
    input = os.path.abspath(input)
    project_name = os.path.basename(input)

    defaults = default_sampleinfo.copy()

    pool = []
    if not os.path.isdir(input):
        click.echo(
            "ERROR - Sequence data folder {} does not exist.".format(project_name)
        )
        ctx.abort()
    elif input != os.getcwd():
        for subfolder in os.listdir(input):
            if os.path.isdir("{}/{}".format(input, subfolder)):
                pool.append(defaults.copy())
                pool[-1]["CG_ID_project"] = project_name
                pool[-1]["CG_ID_sample"] = subfolder
    else:
        project_name = "default_sample_info"
        pool.append(defaults.copy())

    with open("{}/{}.json".format(os.getcwd(), project_name), "w") as output:
        json.dump(pool, output, indent=2)
    click.echo("INFO - Created {}.json in current folder".format(project_name))
    done()


@utils.group()
@click.pass_context
def resync(ctx):
    """Updates internal ST with pubMLST equivalent"""


@resync.command()
@click.option(
    "--type",
    default="list",
    type=click.Choice(["report", "list"]),
    help="Output format",
)
@click.option("--customer", default="all", help="Customer id filter")
@click.option(
    "--skip_update", default=False, help="Skips downloading of references", is_flag=True
)
@click.option(
    "--email",
    default=preset_config["regex"]["mail_recipient"],
    help="Forced e-mail recipient",
)
@click.option("--output", help="Full path to output folder", default="")
@click.pass_context
def review(ctx, type, customer, skip_update, email, output):
    """Generates information about novel ST"""
    # Trace exists by some samples having pubMLST_ST filled in. Make trace function later
    ctx.obj["config"]["regex"]["mail_recipient"] = email
    ext_refs = Referencer(config=ctx.obj["config"], log=ctx.obj["log"])
    if not skip_update:
        ext_refs.update_refs()
        ext_refs.resync()
    click.echo("INFO - Version check done. Generating output")
    if type == "report":
        codemonkey = Reporter(
            config=ctx.obj["config"], log=ctx.obj["log"], output=output
        )
        codemonkey.report(type="st_update", customer=customer)
    elif type == "list":
        ext_refs.resync(type=type)
    done()


@resync.command()
@click.argument("sample_name")
@click.option(
    "--force",
    default=False,
    is_flag=True,
    help="Resolves sample without checking for pubMLST match",
)
@click.pass_context
def overwrite(ctx, sample_name, force):
    """Flags sample as resolved"""
    ext_refs = Referencer(config=ctx.obj["config"], log=ctx.obj["log"])
    ext_refs.resync(type="overwrite", sample=sample_name, ignore=force)
    done()
