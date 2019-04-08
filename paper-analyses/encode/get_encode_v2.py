#!/usr/bin/env python3

# https://www.encodeproject.org/report/?type=Experiment&status=released&assay_title=ChIP-seq&award.project=ENCODE&assembly=hg19&biosample_ontology.classification=cell+line&biosample_ontology.term_name=GM12878
# gm12878_accs.txt
# chose only gm12878 chipseq from ENCODE data search, hg19

import requests, json
import sys

# Force return from the server in JSON format
HEADERS = {'accept': 'application/json'}

# Get accs
accs = [item.strip() for item in open("gm12878_accs.txt", "r").readlines()]

# For each acc:
# get unfiltered bam
# get peaks file

def FileURL(acc, suffix):
    return "https://www.encodeproject.org/files/%s/@@download/%s.%s"%(acc, acc, suffix)

for acc in accs:
    URL="https://www.encodeproject.org/experiments/%s/"%(acc)
    response = requests.get(URL, headers=HEADERS)
    response_json_dict = response.json()
    target = response_json_dict["target"]["title"].split()[0]
    bam_accs = {} # num rep->acc
    peak_accs = {} # num rep->acc
    files = response_json_dict["files"]
    for fdata in files:
        if "assembly" not in fdata or not fdata["assembly"] == "hg19": continue
        if len(fdata["biological_replicates"]) > 1: continue
        try:
            num_rep = fdata["biological_replicates"][0]
        except:
            continue
        if (fdata["file_type"] == "bam" and fdata["output_type"]=="unfiltered alignments"):
            if len(fdata["derived_from"])==3:
                rtype = "Paired"
            elif len(fdata["derived_from"])==2: 
                rtype = "Single"
            else: rtype = "unknown"
            bam_accs[num_rep] =(rtype,fdata["accession"])
        if (fdata["output_type"] in ["peaks", "peaks and background as input for IDR"] and \
           fdata["file_type"] in ["bed narrowPeak", "bed broadPeak"]):
            peak_accs[num_rep] = fdata["accession"]
    for repnum in bam_accs.keys():
        if repnum in peak_accs:
            ba = bam_accs[repnum]
            peak_acc = peak_accs[repnum]
            sys.stdout.write(",".join(["GM12878",target, ba[0], FileURL(ba[1], "bam"), FileURL(peak_acc, "bed.gz")])+"\n")

