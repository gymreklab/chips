#!/usr/bin/env python3

# get ischip params from chips json file

import json
import sys
import numpy as np

model = json.load(open(sys.argv[1], "r"))
pcr = int(sys.argv[2])
highbg = int(sys.argv[3])

k = model["frag"]["k"]
theta = model["frag"]["theta"]
mean = k*theta
var = k*theta**2
sd = np.sqrt(var)

spot = int(model["pulldown"]["s"]*100)
#frac = int(model["pulldown"]["f"]*100+1)

if highbg > 0:
    bg = 4
else: bg = 1

# complains if fragment length var is > 1
# no reads output for pcr 20. setting to 10
sys.stdout.write(" -L %s:%s --pcr %s --ground %s:%s"%(np.log(mean), 0.3, pcr, spot, bg))


