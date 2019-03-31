#!/usr/bin/env python3

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

# GM12878,RCOR1,Single,https://www.encodeproject.org/files/ENCFF489JLD/@@download/ENCFF489JLD.bam,https://www.encodeproject.org/files/ENCFF489ZRF/@@download/ENCFF489ZRF.bed.gz

def FileURL(acc, suffix):
    return "https://www.encodeproject.org/files/%s/@@download/%s.%s"%(acc, acc, suffix)

for acc in accs:
    URL="https://www.encodeproject.org/experiments/%s/"%(acc)
    response = requests.get(URL, headers=HEADERS)
    response_json_dict = response.json()
    target = response_json_dict["target"]["title"].split()[0]
    bam_accs = []
    peak_acc = None
    files = response_json_dict["files"]
    for fdata in files:
        if "assembly" not in fdata or not fdata["assembly"] == "hg19": continue
        num_rep = len(fdata["biological_replicates"])
        if (fdata["file_type"] == "bam" and fdata["output_type"]=="unfiltered alignments"):
            if len(fdata["derived_from"])==3:
                rtype = "Paired"
            elif len(fdata["derived_from"])==2: 
                rtype = "Single"
            else: rtype = "unknown"
            bam_accs.append((rtype,fdata["accession"]))
        if (fdata["output_type"]=="optimal idr thresholded peaks") and \
           fdata["file_type"] in ["bed narrowPeak", "bed broadPeak"]:
            peak_acc = fdata["accession"]
    if peak_acc is not None:
        for ba in bam_accs:
            sys.stdout.write(",".join(["GM12878",target, ba[0], FileURL(ba[1], "bam"), FileURL(peak_acc, "bed.gz")])+"\n")

