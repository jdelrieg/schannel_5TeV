#!/usr/bin/env python
import lz4.frame as lz4f
import pickle 
import joblib
import json
import time
import cloudpickle
import gzip
import os, sys
from optparse import OptionParser

import uproot
import numpy as np
from coffea import hist, processor
from coffea.util import load, save
from coffea.nanoevents import NanoAODSchema

#import tt5TeV
import tchannel5TeV_pruebas
from cafea.modules import samples
from cafea.modules import fileReader

import time

start = time.time()


#import linecache2

#if "#" in linecache2.getline(__file__,17).strip():
#	print("We are doing tchannel \n")
#else:
#	print("We are doing ttbar \n")



if __name__ == '__main__':

  import argparse
  parser = argparse.ArgumentParser(description='You can customize your run')
  parser.add_argument('jsonFiles'           , nargs='?', default=''           , help = 'Json file(s) containing files and metadata')
  parser.add_argument('--prefix', '-r'     , nargs='?', default=''           , help = 'Prefix or redirector to look for the files')
  parser.add_argument('--test','-t'       , action='store_true'  , help = 'To perform a test, run over a few events in a couple of chunks')
  parser.add_argument('--pretend'        , action='store_true'  , help = 'Read json files but, not execute the analysis')
  parser.add_argument('--nworkers','-n'   , default=8  , help = 'Number of workers')
  parser.add_argument('--chunksize','-s'   , default=100000  , help = 'Number of events per chunk')
  parser.add_argument('--nchunks','-c'   , default=None  , help = 'You can choose to run only a number of chunks')
  parser.add_argument('--outname','-o'   , default='plots5TeV', help = 'Name of the output file with histograms')
  parser.add_argument('--outpath','-p'   , default='histos5TeV', help = 'Name of the output directory')
  parser.add_argument('--treename'   , default='Events', help = 'Name of the tree inside the files')
  parser.add_argument('--jobs', '-j',  action='store_true', help = 'send jobs')
  parser.add_argument('--queue', '-q'  , default='batch', help = 'Queue to send jobs')
  parser.add_argument('--exclude', '-x', default=None, help = 'Exclude nodes')
  parser.add_argument('--model', '-m', default=None, help = 'Models')
  parser.add_argument('--split', default=None, help = 'Split json file')
  
  
  args = parser.parse_args()
  jsonFiles        = args.jsonFiles
  prefix           = args.prefix
  dotest           = args.test
  nworkers         = int(args.nworkers)
  chunksize        = int(args.chunksize)
  nchunks          = int(args.nchunks) if not args.nchunks is None else args.nchunks
  outname          = args.outname
  outpath          = args.outpath
  pretend          = args.pretend
  treename         = args.treename
  jobs             = args.jobs
  queue            = args.queue
  exclude          = args.exclude
  pathModels       = args.model

  if dotest:
    nchunks = 1 #2
    chunksize = 200
    nworkers = 1
    print('Running a fast test with %i workers, %i chunks of %i events'%(nworkers, nchunks, chunksize))

  if jobs:
    command = " ".join(sys.argv[:])
    command = command.replace(" --jobs ", " ").replace("-j", " ")
    if not command.startswith('python '): command = 'python ' +  command
    commandJob = 'sbatch -p %s -c %i %s --wrap "%s"'%(queue, nworkers, ("-x " + exclude) if exclude is not None else '', command)
    print(commandJob)
    os.system(commandJob)
    exit()

  import multiprocessing
  nCPUcores = multiprocessing.cpu_count()
  if nworkers > nCPUcores:
    print("WAIT !! You only have %i available cores in this PC!!"%nCPUcores)
    exit()
  elif nworkers == nCPUcores:
    print("WARNING: Your are running with all the cores in your machine!")
  

  ### Load samples from json
  samplesdict = {}
  allInputFiles = []

  def LoadJsonToSampleName(jsonFile, prefix):
    sampleName = jsonFile if not '/' in jsonFile else jsonFile[jsonFile.rfind('/')+1:]
    if sampleName.endswith('.json'): sampleName = sampleName[:-5]
    with open(jsonFile) as jf:
      samplesdict[sampleName] = json.load(jf)
      samplesdict[sampleName]['redirector'] = prefix

  if   isinstance(jsonFiles, str) and ',' in jsonFiles: jsonFiles = jsonFiles.replace(' ', '').split(',')
  elif isinstance(jsonFiles, str)                     : jsonFiles = [jsonFiles]
  for jsonFile in jsonFiles:
    if os.path.isdir(jsonFile):
      if not jsonFile.endswith('/'): jsonFile+='/'
      for f in os.listdir(jsonFile):
        if f.endswith('.json'): allInputFiles.append(jsonFile+f)
    else:
      allInputFiles.append(jsonFile)

  # Read from cfg files
  for f in allInputFiles:
    if not os.path.isfile(f):
      raise Exception('[ERROR] Input file' + f + ' not found!')
    # This input file is a json file, not a cfg
    if f.endswith('.json'): 
      LoadJsonToSampleName(f, prefix)
    # Open cfg files
    else:
      with open(f) as fin:
        print(' >> Reading json from cfg file...')
        lines = fin.readlines()
        for l in lines:
          if '#' in l: l=l[:l.find('#')]
          l = l.replace(' ', '').replace('\n', '')
          if l == '': continue
          if ',' in l:
            l = l.split(',')
            for nl in l:
              if not os.path.isfile(l): prefix = nl
              else: LoadJsonToSampleName(nl, prefix)
          else:
            if not os.path.isfile(l): prefix = l
            else: LoadJsonToSampleName(l, prefix)

  flist = {};
  for sname in samplesdict.keys():
    redirector = samplesdict[sname]['redirector']
    flist[sname] = [(redirector+f) for f in samplesdict[sname]['files']]
    samplesdict[sname]['year'] = samplesdict[sname]['year']
    samplesdict[sname]['xsec'] = float(samplesdict[sname]['xsec'])
    samplesdict[sname]['nEvents'] = int(samplesdict[sname]['nEvents'])
    samplesdict[sname]['nGenEvents'] = int(samplesdict[sname]['nGenEvents'])
    samplesdict[sname]['nSumOfWeights'] = float(samplesdict[sname]['nSumOfWeights'])

    # Print file info
    print('>> '+sname)
    print('   - isData?      : %s'   %('YES' if samplesdict[sname]['isData'] else 'NO'))
    print('   - year         : %s'   %samplesdict[sname]['year'])
    print('   - xsec         : %f'   %samplesdict[sname]['xsec'])
    print('   - histAxisName : %s'   %samplesdict[sname]['histAxisName'])
    print('   - options      : %s'   %samplesdict[sname]['options'])
    print('   - tree         : %s'   %samplesdict[sname]['treeName'])
    print('   - nEvents      : %i'   %samplesdict[sname]['nEvents'])
    print('   - nGenEvents   : %i'   %samplesdict[sname]['nGenEvents'])
    print('   - SumWeights   : %f'   %samplesdict[sname]['nSumOfWeights'])
    print('   - Prefix       : %s'   %samplesdict[sname]['redirector'])
    print('   - nFiles       : %i'   %len(samplesdict[sname]['files']))
    for fname in samplesdict[sname]['files']: print('     %s'%fname)
    
    

  if pretend: 
    print('pretending...')
    exit() 

  import pickle as pkl
  model = None

  if pathModels is None:
    #pathModels = 'analysis/tt5TeV/nn/models/model_04Jul22_07h04m.pkl'
    model2 = '/nfs/fanae/user/andreatf/cafea/cafea/analysis/tt5TeV/models/november/rf3j2b_250_4_allvariablesNewMlb_p2v2.pkl'
    model1='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/training/definitive.pkl'     # ESTE es de momento mi estandar para 1bkg
    #model1='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/pruebillas_2bkg/bibkg.pkl'  # ESTE es de momento mi estandar para 2bkg

    modela='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/training_several/best.pkl' 
    modelb='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/training_several/onlyone.pkl' 
    modelc='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/training_several/allbutone.pkl' 
    modeld='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/training_several/better6.pkl' 
    modele='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/training_several/better6_bibkg.pkl'
    modelf='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/pruebillas_2bkg/bibkg.pkl' 

    model2='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/training_several/better2.pkl' 
    model3='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/training_several/better3.pkl' 
    model4='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/training_several/better4.pkl' 
    model5='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/training_several/better5.pkl' 
    
    modelcheck='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/safecheck/training/safecheck.pkl' 
    modelcheck_old='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/safecheck/training/safecheck_oldsamples.pkl'
    prueba='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/safecheck/training/prueba_1.pkl'

    pathModels = [model1, model2]																#ESTA linea es la que vale cuando me quede con uno (en ppio definitivo era el mejor)
    pathModels = [modelcheck,modelcheck_old,modelc,modeld,modele,modelf,model2,model3,model4,model5,modela,modelb,prueba]        #ESTA linea es para testear todas las diferentes combinaciones de variables
    
    modelloosecuts1='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/lesscuts/training/all_vars.pkl' 
    modelloosecuts2='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/lesscuts/training/varspruned.pkl' 
    #modelloosecuts2='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/lesscuts/training/var_pruned_bibkg.pkl' 
    modelloosecuts3='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/lesscuts/training/best7.pkl' 
   
    pathModels = [modelloosecuts1,modelloosecuts2,modelloosecuts3]     											#ESTA linea para testear lo de tX meeting de relajar cortes y meter esas vars como input al mva

    modeleta1='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/eta3loose/training/all.pkl' 
    modeleta2='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/eta3loose/training/pruned.pkl' 
    #modelloosecuts2='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/lesscuts/training/var_pruned_bibkg.pkl' 
    modeleta3='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/eta3loose/training/best6.pkl' 
    modeleta4='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/eta3loose_morevars/training/b10_besthyper.pkl' 
    
    modeleta5='/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/mva/presented_intX_cuts/training/tXcuts.pkl' #

    pathModels = [modeleta1,modeleta2,modeleta3,modeleta4,modeleta5]  
    
  if pathModels is not None and isinstance(pathModels, str) and not ',' in pathModels:
    with open(pathModels, 'rb') as f:
      model = pickle.load(f)
  elif pathModels is not None and isinstance(pathModels, str):
    pathModels = pathModels.replace(' ', '').split(',')
  if isinstance(pathModels, list):
    model = [pkl.load(open(p, 'rb')) for p in pathModels]
    print('')
  processor_instance = tchannel5TeV_pruebas.AnalysisProcessor(samplesdict, model)  #donde pone tchannel5TeV poner tt5TeV para ttbar
  
  





  # Run the processor and get the output
  tstart = time.time()
  output = processor.run_uproot_job(flist, treename=treename, processor_instance=processor_instance, executor=processor.futures_executor, executor_args={"schema": NanoAODSchema,'workers': nworkers}, chunksize=chunksize, maxchunks=nchunks)
  dt = time.time() - tstart

  nbins = sum(sum(arr.size for arr in h._sumw.values()) for h in output.values() if isinstance(h, hist.Hist))
  nfilled = sum(sum(np.sum(arr > 0) for arr in h._sumw.values()) for h in output.values() if isinstance(h, hist.Hist))
  print("Filled %.0f bins, nonzero bins: %1.1f %%" % (nbins, 100*nfilled/nbins if nbins!=0 else 0,))
  print("Processing time: %1.2f s with %i workers (%.2f s cpu overall)" % (dt, nworkers, dt*nworkers, ))
 
  # Save the output
  if not os.path.isdir(outpath): os.system("mkdir -p %s"%outpath)
  out_pkl_file = os.path.join(outpath,outname+".pkl.gz")
  print("\nSaving output in "+out_pkl_file+"...")
  with gzip.open(out_pkl_file, "wb") as fout:
    cloudpickle.dump(output, fout)
  print("Done!")
  
  end = time.time()
  print('it took',end - start,'seconds')
