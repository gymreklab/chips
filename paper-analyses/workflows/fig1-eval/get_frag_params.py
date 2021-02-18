#!/usr/bin/env python3

# get chipulate fragment length params from chips model file

import json
import sys

model = json.load(open(sys.argv[1], "r"))
k = model["frag"]["k"]
theta = model["frag"]["theta"]
mean = k*theta
var = k*theta**2

sys.stdout.write("--fragment-length %s "%(int(mean)))

