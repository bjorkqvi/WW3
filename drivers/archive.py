#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 17:33:52 2021

@author: janvb
"""
import os
from glob import glob
# !!! This should only be run from the work directory!

if os.path.isfile('wrkdir.txt'):
    raise Exception('!!! This should only be run from the work directory !!!')

import os

with open ("run.index", "r") as f:
    runindex=f.read().replace('\n', '').split()
    
#files=['*.list','*.nml','*.inp','*.ww3','*.nc','*tab*','slurm*.out','*.log']
files=['*']

for file in files:
    ff=glob(file)
    if ff:
        os.system("cp -v %s ../archive/%s/%s/" % (file, runindex[0], runindex[1]))
    
