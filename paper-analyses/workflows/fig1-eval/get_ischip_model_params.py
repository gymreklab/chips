#!/usr/bin/env python3

# get ischip params from chips json file

import json
import sys
import numpy as np

model = json.load(open(sys.argv[1], "r"))
k = model["frag"]["k"]
theta = model["frag"]["theta"]
mean = k*theta
var = k*theta**2
sd = np.sqrt(var)

pcr = 0 # no pcr if the rate is 1
#if model["pcr_rate"] < 1:
#    pcr = 10 # ?

spot = int(model["pulldown"]["s"]*100)
frac = int(model["pulldown"]["f"]*100+1)

# complains if fragment length var is > 1
# no reads output for pcr 20. setting to 10
sys.stdout.write(" -L %s:%s --pcr %s --ground %s:%s"%(np.log(mean), 1, pcr, spot, frac))


