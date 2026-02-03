'''
  Obtain a json file to apply Electron Energy Corrections in coffea from standard Electron Scale/Smear files
  WARNING: as it is, it only works with 5.02 scale corrections... format varies a bit for other periods
'''

import json, os, sys
import argparse

def main():
  parser = argparse.ArgumentParser(description='Create a .json file to use with coffea from ElecES .dat files')
  parser.add_argument('path'  ,  help = 'Path to input file')
  parser.add_argument('--ftype'  , '-t', default='scale'           , help = 'Type = scale/smear')
  parser.add_argument('--outname', '-o', default=None           , help = 'Output name')

  args, unknown = parser.parse_known_args()
  path = args.path
  ftype = args.ftype
  outname = args.outname
  ElectronScaleToJson(path, ftype, outname)
 
  #ElectronScaleToJson('scale.txt', 'scale')

def LoadFile(ipath, ftype = 'scale'):
  dlines = []
  with open(ipath) as f:
    for l in f.readlines():
      if l.startswith('#'): continue
      if '#' in l: l = l[:l.index('#')]
      if '\n' in l: l.replace('\n', '')
      if l == '': continue
      if   ftype == 'scale': dlines.append(ProcessLineScale(l))
      elif ftype == 'smear': dlines.append(ProcessLineSmear(l, var='emean'))
      elif ftype == 'smearrho': dlines.append(ProcessLineSmear(l, var='rho'))
  return dlines

def ProcessFile(ipath, ftype='scale'):
  dlines = LoadFile(ipath, ftype)
  dic = {}; 
  for dB,dR in dlines:
    if len(dB) == 2:
      abseta = dB['absEta']
      r9 = dB['R9']
      if not abseta in dic.keys()        : dic[abseta] = {}
      if not r9     in dic[abseta].keys(): dic[abseta][r9] = {}
      for k in dR: dic[abseta][r9][k] = dR[k]
    elif len(dB) == 3:
      abseta = dB['absEta']
      r9 = dB['R9']
      runNumber = dB['runNumber']
      if not abseta    in dic            .keys(): dic[abseta]                = {}
      if not r9        in dic[abseta]    .keys(): dic[abseta][r9]            = {}
      if not runNumber in dic[abseta][r9].keys(): dic[abseta][r9][runNumber] = {}
      for k in dR: dic[abseta][r9][runNumber][k] = dR[k]
  return dic

def ProcessLineSmear(line, var='emean'):
  ''' Format:
   # category                         Emean  err_Emean  rho     err_rho 
   # absEta_2_2.5-R9_0.940_1.000       6.54  0          0.03625  0.000041
  '''
  dicBins = {}; dicResults = {}
  line = line.replace('\t', ' ')
  while '  ' in line: line = line.replace('  ', ' ')
  cat = line.replace('\n','').split(' ')
  nam = cat[0] 
  emean, erremean, rho, errrho = cat[1:5]
  cats = nam.split('-')
  for c in cats:
    name, xmin, xmax = c.split('_')
    xmin = float(xmin); xmax = float(xmax)
    dicBins[name] = "%s:[%1.2f,%1.2f]"%(name, xmin, xmax)
  dicResults['value'] = float(emean) if var=='emean' else float(rho)
  dicResults['error'] = float(erremean) if var=='emean' else float(errrho)
  return dicBins, dicResults

def ProcessLineScale(line):
  ''' Format: 
    absEta_0_1-R9_0.940_1.000-gainEle_12  runNumber 307042  307063  1.002 0.00059 0.0 0.00086 0.0
  '''
  dicBins = {}; dicResults = {}
  line = line.replace('\t', ' ')
  while '  ' in line: line = line.replace('  ', ' ')
  cat = line.replace('\n','').split(' ')
  nam = cat[0] 
  cats = nam.split('-')
  for c in cats:
    if c.startswith('gain'): continue
    name, xmin, xmax = c.split('_')
    xmin = float(xmin); xmax = float(xmax)
    dicBins[name] = "%s:[%1.2f,%1.2f]"%(name, xmin, xmax)
  cat = cat[1:]
  # Get runNumber
  runMin = 0; runMax = 0; index = 0
  for i in range(len(cat)):
    if cat[i] == 'runNumber': 
      runMin = cat[i+1]
      runMax = cat[i+2]
      index = i
      break
  if runMax != 0: dicBins['runNumber'] = "runNumber:[%s,%s]"%(runMin, runMax)
  cat = cat[index+3:]
  dicResults['value'] = float(cat[0])
  dicResults['error'] = float(cat[1])
  return dicBins, dicResults

def FixDicElecES(dic, var='runNumber'):
  dics = {"value": 1.00, "error": 0.0}
  newDic = dic.copy()
  for d1 in dic:
    for d2 in dic[d1]:
      rmin = []; rmax = []; rundel=[]; newbins=[]; values=[]
      # open intervals in +1 run number
      for runs in dic[d1][d2].keys():
        runsmin, runmax = (runs[runs.index('[')+1:runs.index(']')]).split(',')
        binName = 'runNumber:[%s,%s]'%(str(int(runsmin)-1), str(int(runmax)+1))
        newbins.append(binName)
        values.append(dic[d1][d2][runs])
        rundel.append(runs)
        rmin.append(runsmin)
        rmax.append(runmax)
      for runs in rundel: del dic[d1][d2][runs]
      for b, val in zip(newbins, values): dic[d1][d2][b] = val
      rmin.sort(); rmax.sort()
      for i in range(len(rmin)-1):
        newBinName = 'runNumber:[%s,%s]'%(str(int(rmax[i])+1), str(int(rmin[i+1])-1))
        newDic[d1][d2][newBinName] = dics.copy()
  for k in newDic.keys():
    if ',1.44]' in k:
      rmkey  = k
      addkey = k.replace(',1.44]', ',1.57]')
      val = newDic[k]
      break
  newDic[addkey] = val
  del newDic[rmkey]
  return newDic
        

def ElectronScaleToJson(ifile, ftype='scale', ofile=None):
  if not os.path.isfile(ifile):
    print('Input file does not exist!')
    return
  if ofile is None: ofile = ifile[:ifile.rfind('.')] + '.json'
  dic = ProcessFile(ifile, ftype)
  if ftype=='scale': dic = FixDicElecES(dic)
  dic = {'correction':{ftype:dic}}
  with open(ofile, 'w') as of:
    json.dump(dic, of, indent=2)
    print('>>> New json file: ', ofile)

if __name__ == '__main__':
  main()

