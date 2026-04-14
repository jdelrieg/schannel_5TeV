#!/usr/bin/env python3

import os
from config import *
from multiprocessing import Pool
import ROOT


if var is None: var = "MVAscore_relaxed_b10"

# "path" es la variable que le damos al coso esti.
outpath = path + "combineFiles/"
if not os.path.exists(outpath):
    os.makedirs(outpath)

if outpatho == None:
    raise RuntimeError("FATAL: no output folder set for this cards production!")
pathcomb = path + "/temp_cards/" + outpatho
if not os.path.exists(pathcomb):
    os.makedirs(pathcomb)


levels   = ['2j1b','2j2b','3j2b']
channels = ['e','m']


def getDatacard(task):
    ch, level = task
    print("> Canal:", ch, "Nivel:", level)
    var = "counts" #"MVAscore_relaxed_b10" if level in ['2j1b'] else  "absu0eta"
    outname = "%s_%s_%s.root"%(var, ch, level)
    if not os.path.exists(f"{pathcomb+outname}"):
        command = "python analysis/tt5TeV/SaveRootfile.py -p %s -v %s -l %s -c %s --data"%(path, var, level, ch)  
        if verbose >= 1: print("Running: %s"%(command))
        os.system(command)

        # Move the file to the combine folder
        mvcommand = f"cp {outpath+outname} {pathcomb+outname}"
        if verbose: print("Running: %s"%(mvcommand))
        os.system(mvcommand)

    # Create the datacard
    cardcommand = f"python analysis/tt5TeV/CreateDatacard_v2.py -p {path} --inputFile {pathcomb+outname}"
    if verbose >= 1: print("Running: %s"%(cardcommand))
    os.system(cardcommand)
    if verbose >= 1: print("  ")
    return

tsklist = []
for ch in channels:
    for level in levels:
        tsklist.append( (ch, level) )
print("\n### Iniciando procesamiento")
if nSlots <= 1:
    for tsk in tsklist: getDatacard(tsk)
else:
    pool = Pool(nSlots)
    pool.map(getDatacard, tsklist)
    pool.close()
    pool.join()
print('\n### Datacards created in ', pathcomb, "\n")





