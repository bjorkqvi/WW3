#!/usr/bin/env python3
import os
from glob import glob
import sys

if not len(sys.argv)==3:
    raise Exception('Input: switch and grid, e.g. R33_021 Gslf250v008')

switch=sys.argv[1]
grid=sys.argv[2]

with open ("wrkdir.txt", "r") as f:
    wrkdir=f.read().replace('\n', '')

folder="%s/archive/%s_%s" % (wrkdir, switch,grid)

if not os.path.isdir(folder):
    os.mkdir(folder)
    print ("Creating folder %s..." % folder)
    
dirlist=glob(folder+"/"+"*"+"/")
if dirlist:
    dirlist.sort() # Make sure highest index is last
    oldind=int(dirlist[-1][-4:-1])
    newind=oldind+1
else:
    newind=1

if newind>999:
    raise Exception('Index can only go to 999')

newdir="%s/%03d/" % (folder,newind)
print("Creating new archive for run %s..." % newdir)
os.mkdir(newdir)

# This info is needed for the archiving


with open((('%s/work/run.index') % wrkdir),'w') as f:
    print("Copying archive.py to work directory %s/work and creating the run.index input..." % wrkdir)
    runind="%s_%s %03d\n" % (switch,grid,newind)
    f.write(runind)
    os.system("cp -v archive.py %s/work" % (wrkdir))
    os.system("cp -v convert_tab.sh %s/work" % (wrkdir))

# If this is not the first run with this combination of switches and grids then copy the nml/inp from previous run
if newind>1:
    print("Pulling input files from previous run...")
    os.system("python pull_input.py %s %s %03d" % (switch, grid, oldind))
    
