from datetime import date
from flask import Flask, render_template
from flask_mail import Mail, Message
import logging
from xhtml2pdf import pisa
from io import StringIO, BytesIO

from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import *
from sqlalchemy.sql.expression import case, func

from microSALT.store.db_manipulator import app
from microSALT.store.orm_models import Projects, Samples, Seq_types, Versions

engine = create_engine(app.config['SQLALCHEMY_DATABASE_URI'])
Session = sessionmaker(bind=engine)
session = Session()
mail_ext = Mail()
app.debug = 0
mail_ext.init_app(app)
#Removes server start messages
log = logging.getLogger('werkzeug')
log.setLevel(logging.CRITICAL)

@app.route('/')
def start_page():
    projects = session.query(Projects).all()

    return render_template('start_page.html',
        projects = projects)

@app.route('/microSALT/')
def reroute_page():
    projects = session.query(Projects).all()

    return render_template('start_page.html',
        projects = projects)

@app.route('/microSALT/<project>')
def project_page(project):
    organism_groups = list()
    organism_groups.append("all")
    distinct_organisms = session.query(Samples).filter_by(CG_ID_project=project).distinct()
    for one_guy in distinct_organisms:
        if one_guy.organism not in organism_groups and one_guy.organism is not None:
            organism_groups.append(one_guy.organism)
    organism_groups.sort()
    return render_template('project_page.html',
        organisms = organism_groups,
        project = project) 

@app.route('/microSALT/<project>/<organism_group>')
def report_page(project, organism_group):
    if organism_group == "all":
        sample_info = session.query(Samples).filter(Samples.CG_ID_project==project)
    else:
        sample_info = session.query(Samples).\
        filter(Samples.organism==organism_group, Samples.CG_ID_project==project)
    sample_info = sorted(sample_info, key=lambda sample: int(sample.CG_ID_sample.replace(sample.CG_ID_project, '')[1:]))
    versions = session.query(Versions).all()
    versiondict = dict()
    for version in versions:
      name = version.name[8:]
      versiondict[name] = version.version

    return render_template('report_page.html',
        samples = sample_info,
        date = date.today().isoformat(),
        version = versiondict)
