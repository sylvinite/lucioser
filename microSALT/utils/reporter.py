"""Generates various reports by tapping into Flask and mySQL
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python
import json
import requests
import os
import socket
import sys
import smtplib
import time
import yaml

from datetime import datetime
from shutil import copyfile

from os.path import basename
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication

from multiprocessing import Process

from microSALT import __version__
from microSALT.server.views import app, session, gen_reportdata, gen_collectiondata
from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.store.orm_models import Samples


class Reporter:
    def __init__(
        self, config, log, sampleinfo={}, name="", output="", collection=False
    ):
        self.db_pusher = DB_Manipulator(config, log)
        self.name = name
        self.collection = collection
        if output == "":
            self.output = os.getcwd()
        else:
            self.output = output + "/"
        self.config = config
        self.logger = log
        for k, v in config.items():
            app.config[k] = v
        self.server = Process(target=app.run)
        self.attachments = list()
        self.filedict = dict()
        self.error = False
        self.dt = datetime.now()
        self.now = time.strftime(
            "{}.{}.{}_{}.{}.{}".format(
                self.dt.year,
                self.dt.month,
                self.dt.day,
                self.dt.hour,
                self.dt.minute,
                self.dt.second,
            )
        )

        self.sampleinfo = sampleinfo
        self.sample = None
        if isinstance(self.sampleinfo, list) and len(self.sampleinfo) > 1:
            self.name = self.sampleinfo[0].get("CG_ID_project")
            self.sample = self.sampleinfo[0]
            for entry in self.sampleinfo:
                if entry.get("CG_ID_sample") == self.name:
                    raise Exception(
                        "Mixed projects in samples_info file. Do not know how to proceed"
                    )
        else:
            if isinstance(self.sampleinfo, list):
                self.sampleinfo = self.sampleinfo[0]
            self.name = self.sampleinfo.get("CG_ID_project")
            self.sample = self.sampleinfo

    def create_subfolders(self):
        os.makedirs("{0}/deliverables".format(self.config["folders"]["reports"]), exist_ok=True)
        os.makedirs("{0}/json".format(self.config["folders"]["reports"]), exist_ok=True)
        os.makedirs("{0}/analysis".format(self.config["folders"]["reports"]), exist_ok=True)

    def report(self, type="default", customer="all"):
        self.create_subfolders()
        if type in ["default", "typing", "qc"]:
            # Only typing and qc reports are version controlled
            self.gen_version(self.name)
        if type in ["default", "typing", "qc", "st_update"]:
            self.restart_web()
            if type == "default":
                self.gen_typing()
                self.gen_qc()
                self.gen_json(silent=True)
                self.gen_delivery()
                self.gen_multiqc()
            elif type == "typing":
                self.gen_typing()
            elif type == "qc":
                self.gen_qc()
            elif type == "st_update":
                self.gen_STtracker(customer)
            self.kill_flask()
        elif type in ["json_dump", "motif_overview"]:
            if type == "json_dump":
                self.gen_json()
                self.gen_delivery()
            elif type == "motif_overview":
                self.gen_motif(motif="resistance")
                self.gen_motif(motif="custom")
            elif type == "multiqc":
                self.gen_multiqc()
        else:
            raise Exception("Report function recieved invalid format")
        self.mail()
        #If no output dir is specified; Don't store report locally. Rely on e-mail
        if not self.output == "" or self.output == os.getcwd():
            for k,v in self.filedict.items():
                if v == "": 
                    os.remove(k)
                else:
                    copyfile(k, v)

    def gen_version(self, name):
        self.db_pusher.get_report(name)
        self.db_pusher.set_report(name)

    def gen_STtracker(self, customer="all", silent=False):
        self.name = "Sequence Type Update"
        try:
            r = requests.get(
                "http://127.0.0.1:5000/microSALT/STtracker/{}".format(customer),
                allow_redirects=True,
            )
            outname = "{}/ST_updates_{}.html".format(self.output, self.now)
            outfile = open(outname, "wb")
            outfile.write(r.content.decode("iso-8859-1").encode("utf8"))
            outfile.close()
            self.filedict[outname] = ""
            if not silent:
                self.attachments.append(outname)
        except Exception as e:
            self.logger.error(
                "Flask instance currently occupied. Possible rogue process. Retry command"
            )
            self.error = True

    def gen_typing(self, silent=False):
        try:
            last_version = self.db_pusher.get_report(self.name).version
        except Exception as e:
            self.logger.error("Project {} does not exist".format(self.name))
            self.kill_flask()
            sys.exit(-1)
        try:
            r = requests.get(
                "http://127.0.0.1:5000/microSALT/{}/typing/all".format(self.name),
                allow_redirects=True,
            )
            outfile = "{}_Typing_{}.html".format(
                self.sample.get("Customer_ID_project"), last_version
            )
            local = "{}/{}".format(self.output, outfile)
            output = "{}/analysis/{}".format(self.config["folders"]["reports"], outfile)

            outfile = open(output, "wb")
            outfile.write(r.content.decode("iso-8859-1").encode("utf8"))
            outfile.close()

            if os.path.isfile(output):
                self.filedict[output] = local
                if not silent:
                    self.attachments.append(output)
        except Exception as e:
            self.logger.error(
                "Flask instance currently occupied. Possible rogue process. Retry command"
            )
            self.error = True

    def gen_qc(self, silent=False):
        try:
            last_version = self.db_pusher.get_report(self.name).version
        except Exception as e:
            self.logger.error("Project {} does not exist".format(self.name))
            self.kill_flask()
            sys.exit(-1)
        try:
            q = requests.get(
                "http://127.0.0.1:5000/microSALT/{}/qc".format(self.name),
                allow_redirects=True,
            )
            outfile = "{}_QC_{}.html".format(
                self.sample.get("Customer_ID_project"), last_version
            )
            local = "{}/{}".format(self.output, outfile)
            output = "{}/analysis/{}".format(self.config["folders"]["reports"], outfile)

            outfile = open(output, "wb")
            outfile.write(q.content.decode("iso-8859-1").encode("utf8"))
            outfile.close()

            if os.path.isfile(output):
                self.filedict[output] = local
                if not silent:
                    self.attachments.append(output)
        except Exception as e:
            self.logger.error(
                "Flask instance currently occupied. Possible rogue process. Retry command"
            )
            self.error = True

    def gen_multiqc(self, silent=False):
        try:
            last_version = self.db_pusher.get_report(self.name).version
        except Exception as e:
            self.logger.error("Project {} does not exist".format(self.name))
            sys.exit(-1)
        try:
            outfile = "{}_multiqc_{}.html".format(
                self.sample.get("Customer_ID_project"), last_version
            )
            local = "{}/{}".format(self.output, outfile)
            
            cmd = "multiqc {0} -o {0}/multiqc --no-megaqc-upload".format(self.output)
            proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
            output, error = proc.communicate() 
            cmd = "mv multiqc/multiqc_report.html {0} && rm -r {0}/multiqc".format(local, self.output)
            proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
            output, error = proc.communicate() 

            if os.path.isfile(local):
                self.filedict[output] = local
                if not silent:
                    self.attachments.append(local)
        except Exception as e:
            self.logger.error(
                "Unable to generate multiqc report."
            )
            self.error = True

    def gen_motif(self, motif="resistance", silent=False):
        if motif not in ["resistance", "custom"]:
            self.logger.error("Invalid motif type specified for gen_motif function")
        if self.collection:
            sample_info = gen_collectiondata(self.name)
        else:
            sample_info = gen_reportdata(self.name)
        output = "{}/{}_{}_{}.csv".format(self.output, self.name, motif, self.now)

        # Load motif & gene names into dict
        motifdict = dict()
        for s in sample_info["samples"]:
            if motif == "resistance":
                for r in s.resistances:
                    if (
                        not (r.resistance in motifdict.keys())
                        and r.threshold == "Passed"
                    ):
                        if r.resistance is None:
                            r.resistance = "None"
                        motifdict[r.resistance] = list()
                    if (
                        r.threshold == "Passed"
                        and not r.gene in motifdict[r.resistance]
                    ):
                        motifdict[r.resistance].append(r.gene)
            elif motif == "custom":
                for e in s.custom_targets:
                    if (
                        not (e.virulence in motifdict.keys())
                        and e.threshold == "Passed"
                    ):
                        if e.virulence is None:
                            e.virulence = "None"
                        motifdict[e.virulence] = list()
                    if e.threshold == "Passed" and not e.gene in motifdict[e.virulence]:
                        motifdict[e.virulence].append(e.gene)
        for k, v in motifdict.items():
            motifdict[k] = sorted(v)

        # Top 2 Header
        sepfix = "sep=,"
        topline = "Identity {}% & Span {}%,,,".format(
            self.config["threshold"]["motif_id"], self.config["threshold"]["motif_span"]
        )
        botline = "CG Sample ID,Sample ID,Organism,Sequence Type,Thresholds"
        for k in sorted(motifdict.keys()):
            genes = [""] * len(motifdict[k])
            active_gene = k.replace(",", " &")
            if active_gene == "":
                active_gene = "Uncategorized hits"
            geneholder = ",".join(genes)
            topline += ",,{}{}".format(active_gene, geneholder)
            resnames = ",".join(sorted(motifdict[k]))
            botline += ",,{}".format(resnames)

        try:
            excel = open(output, "w+")
            excel.write("{}\n".format(sepfix))
            excel.write("{}\n".format(topline))
            excel.write("{}\n".format(botline))

            # Create each individual row past the 2nd, per iteration
            for s in sample_info["samples"]:
                rowdict = dict()
                pref = "{},{},{},{},{}".format(
                    s.CG_ID_sample,
                    s.Customer_ID_sample,
                    s.organism,
                    s.ST_status.replace(",", ";"),
                    s.threshold,
                )
                # Load single sample
                if motif == "resistance":
                    for r in s.resistances:
                        if (
                            not (r.resistance in rowdict.keys())
                            and r.threshold == "Passed"
                        ):
                            rowdict[r.resistance] = dict()
                        if (
                            r.threshold == "Passed"
                            and not r.gene in rowdict[r.resistance]
                        ):
                            rowdict[r.resistance][r.gene] = r.identity
                elif motif == "custom":
                    for e in s.custom_targets:
                        if (
                            not (e.virulence in rowdict.keys())
                            and e.threshold == "Passed"
                        ):
                            rowdict[e.virulence] = dict()
                        if (
                            e.threshold == "Passed"
                            and not e.gene in rowdict[e.virulence]
                        ):
                            rowdict[e.virulence][e.gene] = e.identity
                # Compare single sample to all
                hits = ""
                for res in sorted(motifdict.keys()):
                    if res in rowdict.keys():
                        hits += ",1"
                        for gen in sorted(motifdict[res]):
                            hits += ","
                            if gen in rowdict[res].keys():
                                # UPD: Change this to identity of hit
                                hits += "{}".format(rowdict[res][gen])
                            else:
                                hits += "0"
                    else:
                        # Commas eq to res + gen length
                        hits += ",0,0"
                        pad = ["0"] * len(motifdict[res])
                        hits += ",".join(pad)

                excel.write("{}{}\n".format(pref, hits))

            excel.close()
            if os.path.isfile(output):
                self.filedict[output] = ""
                if not silent:
                    self.attachments.append(output)
        except FileNotFoundError as e:
            self.logger.error(
                "Gen_motif unable to produce excel file. Path {} does not exist".format(
                    os.path.basename(output)
                )
            )

    def gen_delivery(self):
        deliv = dict()
        deliv['files'] = list()
        last_version = self.db_pusher.get_report(self.name).version
        output = "{}/deliverables/{}_deliverables.yaml".format(self.config["folders"]["reports"], self.sample.get("Customer_ID_project"))
        local = "{}/{}_deliverables.yaml".format(self.output, self.sample.get("Customer_ID_project"))

        #Project-wide
        #Sampleinfo
        deliv['files'].append({'format':'json','id':self.sample.get("Customer_ID_project"),
                               'path':"{}/sampleinfo.json".format(self.output),
                               'path_index':'~','step':'analysis','tag':'sampleinfo'})
        #QC report
        deliv['files'].append({'format':'html','id':self.sample.get("Customer_ID_project"),
                               'path':"{}/{}_QC_{}.html".format(self.output, self.sample.get("Customer_ID_project"), last_version),
                               'path_index':'~','step':'result_aggregation','tag':'microsalt-qc'})
        #Typing report
        deliv['files'].append({'format':'html','id':self.sample.get("Customer_ID_project"),
                               'path':"{}/{}_Typing_{}.html".format(self.output, self.sample.get("Customer_ID_project"), last_version),
                               'path_index':'~','step':'result_aggregation','tag':'microsalt-type'})
        #Json (vogue) report
        deliv['files'].append({'format':'json','id':self.sample.get("Customer_ID_project"),
                               'path':"{}/{}.json".format(self.output, self.sample.get("CG_ID_project")),
                               'path_index':'~','step':'result_aggregation','tag':'microsalt-json'})
        #Settings dump
        deliv['files'].append({'format':'txt','id':self.sample.get("Customer_ID_project"),
                               'path':"{}/config.log".format(self.output),
                               'path_index':'~','step':'analysis','tag':'runtime-settings'})

        #Sample-wide
        #Single sample
        if self.sampleinfo == self.sample:
            hklist = list()
            hklist.append(self.sampleinfo)
            resultsdir = self.output
        #Project
        else:
            hklist = self.sampleinfo

        for s in hklist:
            if len(hklist) > 1:
                resultsdir = os.path.join(self.output, s["CG_ID_sample"])
            #Contig/Assembly file
            deliv['files'].append({'format':'fasta','id':s["CG_ID_sample"],
                                   'path':"{0}/assembly/{1}_trimmed_contigs.fasta".format(resultsdir, s["CG_ID_sample"]),
                                   'path_index':'~','step':'assembly','tag':'assembly'})
            #Concat trimmed reads forwards
            deliv['files'].append({'format':'fastq','id':s["CG_ID_sample"],
                                   'path':"{0}/trimmed/{1}_trim_front_pair.fastq.gz".format(resultsdir, s["CG_ID_sample"]),
                                   'path_index':'~','step':'concatination','tag':'trimmed-forward-reads'}) 
            #Concat trimmed reads reverse
            deliv['files'].append({'format':'fastq','id':s["CG_ID_sample"],
                                   'path':"{0}/trimmed/{1}_trim_rev_pair.fastq.gz".format(resultsdir, s["CG_ID_sample"]),
                                   'path_index':'~','step':'concatination','tag':'trimmed-reverse-reads'})
            #Concat trimmed reads unpaired
            deliv['files'].append({'format':'fastq','id':s["CG_ID_sample"],
                                   'path':"{0}/trimmed/{1}_trim_unpair.fastq.gz".format(resultsdir, s["CG_ID_sample"]),
                                   'path_index':'~','step':'concatination','tag':'trimmed-unpaired-reads'})            
            #Slurm dump
            deliv['files'].append({'format':'txt','id':s["CG_ID_sample"],
                                   'path':"{0}/slurm_{1}.log".format(resultsdir, s["CG_ID_sample"]),
                                   'path_index':'~','step':'analysis','tag':'logfile'})
            #Quast (assembly) qc report
            deliv['files'].append({'format':'tsv','id':s["CG_ID_sample"],
                                   'path':"{0}/assembly/quast/{1}_report.tsv".format(resultsdir, s["CG_ID_sample"]),
                                   'path_index':'~','step':'assembly','tag':'quast-results'})
            #Alignment (bam, sorted)
            deliv['files'].append({'format':'bam','id':s["CG_ID_sample"],
                                   'path':"{0}/alignment/{1}_{2}.bam_sort".format(resultsdir, s["CG_ID_sample"], s["reference"]),
                                   'path_index':'~','step':'alignment','tag':'reference-alignment-sorted'})
            #Alignment (bam, sorted, deduplicated)
            deliv['files'].append({'format':'bam','id':s["CG_ID_sample"],
                                   'path':"{0}/alignment/{1}_{2}.bam_sort_rmdup".format(resultsdir, s["CG_ID_sample"], s["reference"]),
                                   'path_index':'~','step':'alignment','tag':'reference-alignment-deduplicated'})
            #Picard insert size stats
            deliv['files'].append({'format':'meta','id':s["CG_ID_sample"],
                                   'path':"{0}/alignment/{1}_{2}.stats.ins".format(resultsdir, s["CG_ID_sample"], s["reference"]),
                                   'path_index':'~','step':'insertsize_calc','tag':'picard-insertsize'})


        with open(output, 'w') as delivfile:
            documents = yaml.dump(deliv, delivfile)

        with open(output, 'r') as delivfile:
            postfix = delivfile.read()
        postfix = postfix.replace("'~'", "~")

        with open(output, 'w') as delivfile:
            delivfile.write(postfix)

        if os.path.isfile(output):
            self.filedict[output] = local


    def gen_json(self, silent=False):
        report = dict()
        local = "{}/{}.json".format(self.output, self.name)
        output = "{}/json/{}.json".format(self.config["folders"]["reports"], self.name)

        sample_info = gen_reportdata(self.name)
        analyses = [
            "blast_pubmlst",
            "quast_assembly",
            "blast_resfinder_resistence",
            "picard_markduplicate",
            "microsalt_samtools_stats",
        ]
        for s in sample_info["samples"]:
            t = dict()

            # Since some apps are too basic to filter irrelevant non-standard values..
            t["ST_status"] = (
                "" if s.ST_status is None or s.ST_status != str(s.ST) else s.ST_status
            )
            t["threshold"] = (
                ""
                if s.threshold is None or s.threshold not in ["Passed", "Failed"]
                else s.threshold
            )
            t["genome_length"] = (
                ""
                if s.genome_length is None or s.genome_length < 1
                else s.genome_length
            )
            t["gc_percentage"] = (
                ""
                if s.gc_percentage is None or s.gc_percentage < 0.1
                else str(s.gc_percentage)
            )
            t["n50"] = "" if s.n50 is None or s.n50 < 1 else s.n50
            t["contigs"] = "" if s.contigs is None or s.contigs < 1 else s.contigs
            t["insert_size"] = (
                "" if s.insert_size is None or s.insert_size < 1 else s.insert_size
            )
            t["duplication_rate"] = (
                ""
                if s.duplication_rate is None or s.duplication_rate < 0.1
                else s.duplication_rate
            )
            t["total_reads"] = (
                "" if s.total_reads is None or s.total_reads < 1 else s.total_reads
            )
            t["mapped_rate"] = (
                "" if s.mapped_rate is None or s.mapped_rate < 0.1 else s.mapped_rate
            )
            t["average_coverage"] = (
                ""
                if s.average_coverage is None or s.average_coverage < 0.1
                else s.average_coverage
            )
            t["coverage_10x"] = (
                "" if s.coverage_10x is None or s.coverage_10x < 0.1 else s.coverage_10x
            )
            t["coverage_30x"] = (
                "" if s.coverage_30x is None or s.coverage_30x < 0.1 else s.coverage_30x
            )
            t["coverage_50x"] = (
                "" if s.coverage_50x is None or s.coverage_50x < 0.1 else s.coverage_50x
            )
            t["coverage_100x"] = (
                ""
                if s.coverage_100x is None or s.coverage_100x < 0.1
                else s.coverage_100x
            )

            report[s.CG_ID_sample] = dict()
            for a in analyses:
                if a == "blast_resfinder_resistence":
                    report[s.CG_ID_sample][a] = list()
                else:
                    report[s.CG_ID_sample][a] = dict()

            report[s.CG_ID_sample]["blast_pubmlst"] = {
                "sequence_type": t["ST_status"],
                "thresholds": t["threshold"],
            }
            report[s.CG_ID_sample]["quast_assembly"] = {
                "estimated_genome_length": t["genome_length"],
                "gc_percentage": t["gc_percentage"],
                "n50": t["n50"],
                "necessary_contigs": t["contigs"],
            }
            report[s.CG_ID_sample]["picard_markduplicate"] = {
                "insert_size": t["insert_size"],
                "duplication_rate": t["duplication_rate"],
            }
            report[s.CG_ID_sample]["microsalt_samtools_stats"] = {
                "total_reads": t["total_reads"],
                "mapped_rate": t["mapped_rate"],
                "average_coverage": t["average_coverage"],
                "coverage_10x": t["coverage_10x"],
                "coverage_30x": t["coverage_30x"],
                "coverage_50x": t["coverage_50x"],
                "coverage_100x": t["coverage_100x"],
            }

            for r in s.resistances:
                if (
                    not (r.gene in report[s.CG_ID_sample]["blast_resfinder_resistence"])
                    and r.threshold == "Passed"
                ):
                    report[s.CG_ID_sample]["blast_resfinder_resistence"].append(r.gene)

        # json.dumps(report) #Dumps the json directly
        try:
            with open(output, "w") as outfile:
                json.dump(report, outfile)

            if os.path.isfile(output):
                self.filedict[output] = local
                if not silent:
                    self.attachments.append(output)
        except FileNotFoundError as e:
            self.logger.error(
                "Gen_json unable to produce json file. Path {} does not exist".format(
                    os.path.basename(output)
                )
            )

    def mail(self):
        msg = MIMEMultipart()
        if not self.error and self.attachments:
            msg["Subject"] = "{} ({}) Reports".format(
                self.name, self.attachments[0].split("_")[0]
            )
        else:
            msg["Subject"] = "{} Failed Generating Report".format(self.name)

        sender = socket.gethostname()
        sender_fixed = "{}.com".format(os.path.splitext(sender)[0])
        msg["From"] = sender_fixed

        msg["To"] = self.config["regex"]["mail_recipient"]

        if not self.error:
            for file in self.attachments:
                part = MIMEApplication(open(file).read())
                part.add_header(
                    "Content-Disposition",
                    'attachment; filename="%s"' % os.path.basename(file),
                )
                msg.attach(part)

        s = smtplib.SMTP("localhost")
        s.connect()
        s.sendmail(msg["From"], msg["To"], msg.as_string())
        s.quit()
        self.logger.info(
            "Mail containing report sent to {} from {}".format(msg["To"], msg["From"])
        )

    def start_web(self):
        self.server.start()
        self.logger.info("Started webserver on http://127.0.0.1:5000/")
        # Hinders requests before server goes up
        time.sleep(0.15)

    def kill_flask(self):
        self.server.terminate()
        self.server.join()
        self.logger.info("Closed webserver on http://127.0.0.1:5000/")

    def restart_web(self):
        try:
          self.kill_flask()
        except Exception as e:
          pass
        self.start_web()
