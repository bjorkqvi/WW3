#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 17:04:29 2021

@author: janvb
"""

# E.g. python pull_input.py R33_001 Gslf250v001 003

import os
import sys

switch=sys.argv[1]
grid=sys.argv[2]
runind=sys.argv[3]

with open ("wrkdir.txt", "r") as f:
    wrkdir=f.read().replace('\n', '')

folder="%s/archive/%s_%s/%s" % (wrkdir, switch,grid,runind)

os.system("cp -v %s/restart.ww3 %s/work/" % (folder, wrkdir))
os.system("cp -v %s/wind.ww3 %s/work/" % (folder, wrkdir))
os.system("cp -v %s/mod_def.ww3 %s/work/" % (folder, wrkdir))
#os.system("cp -v %s/mod_def.ww3 %s/work/" % (folder, wrkdir))

