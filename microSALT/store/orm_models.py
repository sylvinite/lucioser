"""Samples table definition
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import *
from sqlalchemy.orm import relationship

from microSALT import app

db = SQLAlchemy(app)

class Samples(db.Model):
  __tablename__ = 'samples'
  seq_types = relationship("Seq_types", back_populates="samples")
  projects = relationship('Projects', back_populates='samples')
  resistances = relationship("Resistances", back_populates="samples")
  #steps = relationship("Steps", back_populates="samples")
  expacs = relationship("Expacs", back_populates="samples")
  core_seq_types = relationship("Core_seq_types", back_populates="samples")

  CG_ID_sample = db.Column(db.String(15), primary_key=True, nullable=False)
  CG_ID_project = db.Column(db.String(15), ForeignKey('projects.CG_ID_project'))
  Customer_ID_sample = db.Column(db.String(40))
  organism = db.Column(db.String(30))
  ST = db.Column(db.SmallInteger, default=-1)
  pubmlst_ST = db.Column(db.SmallInteger, default=-1)
  cgST = db.Column(db.SmallInteger, default=-1)
  date_analysis = db.Column(db.DateTime)
  genome_length = db.Column(db.Integer, default=-1)
  gc_percentage = db.Column(db.Float(3,2), default = 0.0)
  n50 = db.Column(db.Integer, default=-1)
  contigs = db.Column(db.Integer, default=-1)
  priority = db.Column(db.String(20))

  total_reads = db.Column(db.Integer) #Fetch from bcl2fastq
  insert_size = db.Column(db.Integer)
  duplication_rate = db.Column(db.Float)
  mapped_rate = db.Column(db.Float)
  coverage_10x = db.Column(db.Float)
  coverage_30x = db.Column(db.Float)
  coverage_50x = db.Column(db.Float)
  coverage_100x = db.Column(db.Float)
  average_coverage = db.Column(db.Float)
  reference_genome = db.Column(db.String(32))

  application_tag = db.Column(db.String(15))
  date_arrival = db.Column(db.DateTime)
  date_analysis = db.Column(db.DateTime)
  date_sequencing = db.Column(db.DateTime)
  date_libprep = db.Column(db.DateTime)
  method_sequencing = db.Column(db.String(15))
  method_libprep = db.Column(db.String(15))

class Seq_types(db.Model):
  __tablename__ = 'seq_types'
  samples = relationship('Samples', back_populates='seq_types')

  CG_ID_sample = db.Column(db.String(15), ForeignKey('samples.CG_ID_sample'), primary_key=True)
  loci = db.Column(db.String(10), primary_key=True)
  allele = db.Column(db.SmallInteger)
  contig_name = db.Column(db.String(20), primary_key=True)
  contig_length = db.Column(db.Integer)
  contig_coverage = db.Column(db.Float(6,2))
  identity = db.Column(db.Float(3,2), default= 0.0)
  span = db.Column(db.Float(3,2), default= 0.0)
  evalue = db.Column(db.String(10))
  bitscore = db.Column(db.SmallInteger)
  subject_length = db.Column(db.Integer)
  st_predictor = db.Column(db.Boolean, default = 0)
  contig_start = db.Column(db.Integer)
  contig_end = db.Column(db.Integer)

class Core_seq_types(db.Model):
  __tablename__ = 'core_seq_types'
  samples = relationship('Samples', back_populates='core_seq_types')

  CG_ID_sample = db.Column(db.String(15), ForeignKey('samples.CG_ID_sample'), primary_key=True)
  loci = db.Column(db.String(10), primary_key=True)
  allele = db.Column(db.SmallInteger)
  contig_name = db.Column(db.String(20), primary_key=True)
  contig_length = db.Column(db.Integer)
  contig_coverage = db.Column(db.Float(6,2))
  identity = db.Column(db.Float(3,2), default= 0.0)
  span = db.Column(db.Float(3,2), default= 0.0)
  evalue = db.Column(db.String(10))
  bitscore = db.Column(db.SmallInteger)
  subject_length = db.Column(db.Integer)
  contig_start = db.Column(db.Integer)
  contig_end = db.Column(db.Integer)

class Resistances(db.Model):
  __tablename__ = 'resistances'
  samples = relationship('Samples', back_populates='resistances')

  CG_ID_sample = db.Column(db.String(15), ForeignKey('samples.CG_ID_sample'), primary_key=True)
  gene = db.Column(db.String(50), primary_key=True)
  instance = db.Column(db.String(30), primary_key=True)
  contig_name = db.Column(db.String(20), primary_key=True)
  contig_length = db.Column(db.Integer)
  contig_coverage = db.Column(db.Float(6,2))
  identity = db.Column(db.Float(3,2), default= 0.0)
  span = db.Column(db.Float(3,2), default= 0.0)
  evalue = db.Column(db.String(10))
  bitscore = db.Column(db.SmallInteger)
  subject_length = db.Column(db.Integer)
  reference = db.Column(db.String(40))
  resistance = db.Column(db.String(120))
  contig_start = db.Column(db.Integer)
  contig_end = db.Column(db.Integer)

class Expacs(db.Model):
  __tablename__ = 'expacs'
  samples = relationship('Samples', back_populates='expacs')

  CG_ID_sample = db.Column(db.String(15), ForeignKey('samples.CG_ID_sample'), primary_key=True)
  gene = db.Column(db.String(50), primary_key=True)
  instance = db.Column(db.String(30), primary_key=True)
  contig_name = db.Column(db.String(20), primary_key=True)
  contig_length = db.Column(db.Integer)
  contig_coverage = db.Column(db.Float(6,2))
  identity = db.Column(db.Float(3,2), default= 0.0)
  span = db.Column(db.Float(3,2), default= 0.0)
  evalue = db.Column(db.String(10))
  bitscore = db.Column(db.SmallInteger)
  subject_length = db.Column(db.Integer)
  reference = db.Column(db.String(40))
  virulence = db.Column(db.String(120))
  contig_start = db.Column(db.Integer)
  contig_end = db.Column(db.Integer)

class Projects(db.Model):
  __tablename__ = 'projects'
  samples = relationship('Samples', back_populates='projects')
  reports = relationship('Reports', back_populates='projects')

  CG_ID_project = db.Column(db.String(15), primary_key=True, nullable=False)
  Customer_ID_project = db.Column(db.String(15))
  date_ordered = db.Column(db.DateTime)
  Customer_ID = db.Column(db.String(15))

class Versions(db.Model):
  __tablename__ = 'versions'

  name = db.Column(db.String(45), primary_key=True, nullable=False)
  version = db.Column(db.String(10))

#Keeps and aggregate step string, makes a new version whenever one is not found
class Reports(db.Model):
  __tablename__ = 'reports'
  projects = relationship('Projects', back_populates='reports')

  CG_ID_project = db.Column(db.String(15), ForeignKey('projects.CG_ID_project'), primary_key=True)
  steps_aggregate = db.Column(db.String(100))
  date = db.Column(db.DateTime)
  version = db.Column(db.Integer, default=1, primary_key=True)

class Collections(db.Model):
  __tablename__ = 'collections'

  ID_collection = db.Column(db.String(15), primary_key=True)
  CG_ID_sample = db.Column(db.String(15), primary_key=True)

#Multi-date support for libprep/sequencing/analysis
#class Steps(db.Model):
#  __tablename__ = 'steps'
#  samples = relationship("Samples", back_populates="steps")
#
#  CG_ID_sample = db.Column(db.String(15), ForeignKey('samples.CG_ID_sample'), primary_key=True)
#  step = db.Column(db.String(40), primary_key=True)
#  method = db.Column(db.String(40), primary_key=True)
#  date = db.Column(db.DateTime)
