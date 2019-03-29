#!/usr/bin/env python3
"""
Print chipmunk model in tab format
"""

import json
import os
import sys

try:
    mfile = sys.argv[1]
except:
    sys.stderr.write(__doc__)

m = json.load(open(mfile, "r"))
k = "%0.3f"%float(m["frag"]["k"])
th = "%0.3f"%float(m["frag"]["theta"])
frac = "%0.3f"%float(m["pulldown"]["f"])
spot = "%0.3f"%float(m["pulldown"]["s"])
pcr = "%0.3f"%float(m["pcr_rate"])
meta = os.path.basename(mfile).split(".")[0].split("_")
celltype = meta[0]
fact = meta[1]
encbam = meta[2]
encbed = meta[3]

sys.stdout.write(",".join([celltype, fact, encbam, encbed, k, th, frac, spot, pcr, "s3://chipmunk-encode-models/%s"%os.path.basename(mfile)])+"\n")
