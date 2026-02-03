import os
from config import *

if var is None: var = "medianDRjj"

outpath = path + "combineFiles_check/"

if not os.path.exists(outpath):
    os.makedirs(outpath)
#datatoday='12apr23'
pathcomb = "../../CMSSW_11_3_4/src/combinescripts/tt5TeV2/"+datatoday + "/" #+var+ "/"
if not os.path.exists(pathcomb):
    os.makedirs(pathcomb)

levels = ['3j1b', '3j2b', '4j1b', '4j2b', 'g5j1b', 'g5j2b']
channels = ['e', 'm']
iteration = 0
total = len(channels)*len(levels)

for ch in channels:
    for level in levels:
        #print("\r[{:<100}] {:.2f} % ".format('#' * int(float(iteration)/total*100), float(iteration)/total*100),end='')
        progress_ratio = float(iteration) / total
        progress_percentage = progress_ratio * 100
        progress_bar = '#' * int(progress_ratio * 100)
        formatted_output = "\r[{:<100}] {:.2f}% ".format(progress_bar, progress_percentage)
        print(formatted_output)
        
        var = "medianDRjj" if level not in ['3j1b'] else "MVAscore" #quitar descomentar 
        outname = "%s_%s_%s.root"%(var, ch, level)
        if not os.path.exists("{}{}".format(pathcomb, outname)):#if not os.path.exists(f"{pathcomb+outname}"):
       
          command = "python analysis/tt5TeV/SaveRootfile.py -p %s -v %s -l %s -c %s --data"%(path, var, level, ch)
          if verbose >= 1: print("Running: %s"%(command))
          os.system(command)

          # Move the file to the combine folder
          mvcommand = "cp {}{} {}{}".format(outpath, outname, pathcomb, outname)#mvcommand = f"cp {outpath+outname} {pathcomb+outname}"
          if verbose: print("Running: %s"%(mvcommand))
          os.system(mvcommand)

        # Create the datacard
        cardcommand = "python analysis/tt5TeV/CreateDatacard.py -p {} --inputFile {}{}".format(path, pathcomb, outname)#cardcommand = f"python analysis/tt5TeV/CreateDatacard.py -p {path} --inputFile {pathcomb+outname}"
        if verbose >= 1: print("Running: %s"%(cardcommand))
        os.system(cardcommand)
        if verbose >= 1: print("  ")
        iteration+= 1
print("\r[{:<100}] {:.2f} % ".format('#' * int(float(iteration)/total*100), float(iteration)/total*100))
print('Datacards created in ', pathcomb)
