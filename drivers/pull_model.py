#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 17:19:20 2021

@author: janvb
"""
import os
import sys

with open ("wrkdir.txt", "r") as f:
    wrkdir=f.read().replace('\n', '')
if len(sys.argv) > 1:

	
	with open (f"{wrkdir}/work/run.index", "r") as f:
	    runindex=f.read().replace('\n', '').split()

	if int(runindex[1])>1:
		files=['ww3_grid','ww3_prnc','ww3_strt','ww3_shel','ww3_ounf','ww3_outp','ww3_bounc','ww3_aux.job','ww3_shel.job']
		print(f"Pulling model files from /archive/{runindex[0]}/{int(runindex[1])-1:03d}/")
		for file in files:
		    os.system("cp -v %s/archive/%s/%03d/%s %s/work/" % (wrkdir, runindex[0], int(runindex[1])-1, file, wrkdir))


else:
	files=['grid','prnc','strt','shel','ounf','outp','ounp','bounc']

	for file in files:
	    os.system("cp -v ../model/exe/ww3_%s %s/work/" % (file, wrkdir))

	jobfiles=['ww3_aux','ww3_shel']

	for jobfile in jobfiles:
	    os.system("cp -v ../job/%s.job %s/work/" % (jobfile, wrkdir))
