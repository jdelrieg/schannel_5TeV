# python copyFromJson.py "cafea/json/signal_samples/private_UL/UL18_ttHJet_b1.json, cafea/json/signal_samples/private_UL/UL18_ttHJet_b2.json, cafea/json/signal_samples/private_UL/UL18_ttHJet_b3.json" -o /pool/phedex/userstorage/juanr/ULnanoAOD/noskim/nanoAODv9/2018/EFTsignals/ttH/

# python copyFromJson.py /TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM --DAS -o /pool/phedex/userstorage/juanr/ULnanoAOD/noskim/nanoAODv9/2018/BLAH/ -p root://cms-xrd-global.cern.ch/

import json
import argparse
import os, sys
from cafea.modules.DASsearch import GetFilesFromDatasets

parser = argparse.ArgumentParser(description='Copy files in json')
parser.add_argument('path'              , default=''           , help = 'Path to json')

parser.add_argument('--DAS'     , action='store_true'           , help = 'Get path from dataset in DAS')
parser.add_argument('--prefix','-p'     , default='root://deepthought.crc.nd.edu/'           , help = 'Prefix to add to the path (e.g. redirector)')
parser.add_argument('--output','-o'     , default='.'          , help = 'Output folder')
parser.add_argument('-n'                , default=-1           , help = 'Number of files to copy')
parser.add_argument('--force', '-f'     , action="store_true" , help = 'Force')
parser.add_argument('--fix',          action="store_true" , help = 'Compare size of input and output files and copy if size difference > threshold (1%)')
parser.add_argument('--pretend'     , action="store_true" , help = 'Print the command but not execute')

args, unknown = parser.parse_known_args()


path = args.path
pref = args.prefix
out  = args.output
isDAS = args.DAS
n    = int(args.n)
force = args.force
fix = args.fix
pretend = args.pretend
files = []

def GetFilesInJson(jf):
  with open(jf) as jsfile:
    js = json.load(jsfile)
  return js['files']

def copyFile(fname, prefix='', output='.'):
  command = 'xrdcp %s%s %s'%(prefix, fname, output)
  if not os.path.isdir(output): os.system('mkdir %s'%output)
  outfile = output+fname[fname.rfind('/')+1:]
  redo = False
  if fix:
    f1_size = GetFileSize(fname, prefix)
    f2_size = GetFileSize(outfile)
    if (f1_size-f2_size)/f1_size > 0.01:
      redo = True
  if os.path.isfile(outfile):
    if fix:
      if redo:
        print("Fixing file: ", outfile, "...")
        print(">> Remote size [%i], local size [%i] (difference: %i %s)" %(f1_size, f2_size, ((f1_size-f2_size)/f1_size*100), '%'))
      else: return
    elif not force: 
      print('File already exists! Skipping...')
      return
    else:
      print('File already exists! Forcing!')
  print('Executing: "%s"...'%command)
  if pretend: return
  os.system(command)

def GetFileSize(fname, prefix=''):
  if prefix != '': command = "xrdfs %s stat %s"%(prefix, fname)
  else: command = 'stat %s'%fname
  out = os.popen(command).read()
  n_size = out.find('Size:')
  size = out[n_size+5:].split('\n')[0].split()[0]
  return int(size)

#GetFileInfo('/store/user/kmohrman/FullProduction/FullR2/UL18/Round1/Batch2/naodOnly_step/v1/nAOD_step_tllq4fNoSchanWNoHiggs0p_all22WCsStartPtCheckV2dim6TopMay20GST_run0/NAOD-00000_7156.root', 'root://deepthought.crc.nd.edu/');
#exit();

if ',' in path: path = path.replace(' ', '').split(',')
else          : path = [path]

for jsonfile in path:
  if not isDAS: files += GetFilesInJson(jsonfile)
  else:         files += GetFilesFromDatasets(jsonfile)

for f in files:
  copyFile(f, pref, out)

