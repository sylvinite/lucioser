#!/usr/bin/env python
from microSALT import __version__
from setuptools import setup, find_packages

version = __version__

try:
    with open("requirements.txt", "r") as f:
        install_requires = [x.strip() for x in f.readlines()]
except IOError:
    install_requires = []

setup(
    name="microSALT",
    version=version,
    long_description=__doc__,
    url="https://github.com/Clinical-Genomics/microSALT",
    author="Isak Sylvin",
    author_email='isak.sylvin@scilifelab.se',
    install_requires=install_requires,
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    data_files=[('logos', ['artwork/microsalt.jpg', 'artwork/swedac.jpg']),
                ('testdata', ['tests/testdata/sampleinfo_samples.json','tests/testdata/sampleinfo_mlst.json','tests/testdata/sampleinfo_projects.json','tests/testdata/sampleinfo_reports.json','tests/testdata/sampleinfo_resistance.json','tests/testdata/sampleinfo_reports.json','tests/testdata/sampleinfo_custom.json','tests/testdata/quast_results.tsv', 'tests/testdata/blast_single_resistance.txt','tests/testdata/blast_single_loci.txt','tests/testdata/alignment.stats.ref','tests/testdata/alignment.stats.raw','tests/testdata/alignment.stats.map','tests/testdata/alignment.stats.ins','tests/testdata/alignment.stats.dup','tests/testdata/alignment.stats.cov','configExample.json']),
                ('testproject/AAA1234_2000.1.2_3.4.5/AAA1234A1/', ['tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A1/dummysequences_1.fastq.gz', 'tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A1/dummysequences_2.fastq.gz']), 
                ('testproject/AAA1234_2000.1.2_3.4.5/AAA1234A2/', ['tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A2/dummysequences_1.fastq.gz', 'tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A2/dummysequences_2.fastq.gz']),
                ('testproject/AAA1234_2000.1.2_3.4.5/AAA1234A3/', ['tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A3/dummysequences_1.fastq.gz', 'tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A3/dummysequences_2.fastq.gz']),
                ('testproject/AAA1234_2000.1.2_3.4.5/AAA1234A4/', ['tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A4/dummysequences_1.fastq.gz', 'tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A4/dummysequences_2.fastq.gz'])],
    entry_points={
        'console_scripts': ['microSALT=microSALT.cli:root'],
    },
)

