# Create datacard from simple "combine" rootfile, using uproot
# Usage: 
#   python cafea/modules/CreateDatacardFromRootfile.py NBtags.root -s tt -b "tW, WJets, QCD, DY" -l 0.015 --bkgUnc "0.2, 0.2, 0.2, 0.2" -u "lepSF, trigSF, btagSF, FSR, ISR" -o test.dat
import uproot, os
import numpy as np
from coffea import hist, processor
from cafea.modules.FixNegValuesRoot import FixNegValuesRoot


class Datacard:
  def __init__(self, fname, signal='', bkgList=[], lumiUnc=0.0, bkgUnc=[], systList=[], nSpaces=10, outname=None, rmNegBins=True, verbose=1):
    if rmNegBins:
      FixNegValuesRoot(fname)
    self.LoadFile(fname)
    self.SetSignal(signal)
    self.SetBkg(bkgList)
    self.SetLumiUnc(lumiUnc)
    self.SetBkgUnc(bkgUnc)
    self.SetSystematics(systList)
    self.nBins = 1
    self.chName = 'ch'
    self.datacards = ''
    self.nSpaces = nSpaces
    self.extraUnc = {}
    self.outname = 'temp.txt'
    self.outpath = './'
    self.SetOutName(outname)
    self.SetVerbose(verbose)

  def LoadFile(self, fname):
    self.fname = fname
    self.f = uproot.open(fname)

  def Close(self):
    self.f.close()

  def SetVerbose(self, verbose):
    self.verbose = verbose

  def SetOutPath(self, outpath):
    if not outpath.endswith('/'): outpath += '/'
    if not os.path.isdir(outpath): os.system('mkdir -p %s'%outpath)
    self.outpath = outpath

  def SetOutName(self, outname):
    if outname is None: return
    self.outname = outname

  def GetOutName(self):
    return self.outpath + self.outname
  
  def SetSignal(self, signal):
    self.signal = signal

  def SetBkg(self, bkg):
    if isinstance(bkg, str) and ',' in bkg: bkg = bkg.replace(' ', '').split(',')
    self.bkg = bkg

  def SetLumiUnc(self, lunc):
    self.lumiUnc = lunc

  def SetSystematics(self, syst):
    if isinstance(syst, str) and ',' in syst: syst = syst.replace(' ', '').split(',')
    self.syst = syst

  def SetBkgUnc(self, bkgUnc):
    if isinstance(bkgUnc, dict): self.bkgUnc = bkgUnc
    else:
      self.bkgUnc = {}
      for i in range(len(bkgUnc)):
        if i < len(self.bkg): self.bkgUnc[self.bkg[i]] = bkgUnc[i]
    for b in self.bkg:
      if not b in self.bkgUnc.keys(): self.bkgUnc[b] = 1

  def AddExtraUnc(self, name, val, pr):
    if   isinstance(pr, str) and ',' in pr: pr = pr.replace(' ', '').split(',')
    elif isinstance(pr, str): pr = [pr]
    self.extraUnc[name] = {'process':pr, 'value':val}

  def nProcess(self):
    return len(self.bkg)+1

  def ProcessList(self):
    return [self.signal] + self.bkg

  def FixStringSpace(self, s):
    while len(s) < self.nSpaces: s+= ' ' 
    return s+' '

  def FixListSpaces(self, l):
    t = ''
    for e in l: t+= self.FixStringSpace(str(e))
    return t

  #########################################################################################
  def AddLine(self, l):
    self.datacard += l + '\n'

  def AddSep(self):
    self.AddLine('------------------------------------')

  def AddTxtHeader(self):
    self.AddLine('imax %i number of bins'%self.nBins)
    self.AddLine('jmax %i processes minuns 1'%len(self.bkg))
    self.AddLine('kmax * number of nuisance parameters')
    self.AddSep()
    self.AddLine('shapes * %s %s $PROCESS $PROCESS_$SYSTEMATIC'%(self.chName, self.fname))
    self.AddSep()
    
  def AddTxtObs(self):
    self.AddLine(self.FixListSpaces(['bin', self.chName]))
    self.AddLine(self.FixListSpaces(['observation', '%i'%int(self.Yield('data_obs'))]))
    self.AddSep()

  def AddTxtRates(self):
    self.AddLine(self.FixListSpaces(['bin',     ''] + [self.chName]*self.nProcess()))
    self.AddLine(self.FixListSpaces(['process', ''] + self.ProcessList()))
    self.AddLine(self.FixListSpaces(['process', ''] + list(range(self.nProcess()))))
    self.AddLine(self.FixListSpaces(['rate',    ''] + self.GetProcessRates()))
    self.AddSep()

  def AddLumiUncTxt(self):
    datacard = self.AddLine(self.FixListSpaces(["Lumi", "lnN"] + ([1+self.lumiUnc]*self.nProcess())))

  def AddBkgNormTxt(self):
    for pr in self.bkg:
      self.AddLine(self.FixListSpaces([pr, 'lnN']) + self.GetNormUncProcess(pr))

  def AddSystUncTxt(self):
    for syst in self.syst:
      self.AddLine( self.FixListSpaces([syst, 'shape']) + self.GetSystUncLine(syst))

  def AddExtraUncTxt(self):
    for syst in self.extraUnc.keys():
      val = self.extraUnc[syst]['value']
      process = self.extraUnc[syst]['process']
      prstring = []
      for ipr in self.ProcessList():
        if ipr in process: prstring.append(1+float(val))
        else:              prstring.append('-')
      self.AddLine( self.FixListSpaces([syst, 'lnN'] + prstring) )
 
  def AddAutoMCStatTxt(self):
    self.AddLine('%s autoMCStats 0 0 1'%self.chName)

  def GetDatacard(self):
    self.datacard = ''
    self.AddTxtHeader()
    self.AddTxtObs()
    self.AddTxtRates()
    self.AddLumiUncTxt()
    self.AddBkgNormTxt()
    self.AddSystUncTxt()
    self.AddExtraUncTxt()
    self.AddAutoMCStatTxt()
    return self.datacard

  def Save(self, fout=None):
    self.SetOutName(fout)
    fout = self.GetOutName()
    if os.path.isfile(fout): os.system('mv %s %s_old'%(fout, fout))
    with open(fout, 'w') as f:
      f.write(self.GetDatacard())
    self.Close()
    if self.verbose: print('Created datacard: ', fout)

  #########################################################################################

  def Yield(self, prName):
    vals, bins = self.f[prName].to_numpy(flow=False) # xxx
    return sum(vals)

  def GetProcessRates(self):
    return ['%1.5f'%(self.Yield(pr)) for pr in self.ProcessList()]

  def GetNormUncProcess(self, pr):
    prstring = []
    for ipr in self.ProcessList():
      if ipr == pr: prstring.append(1+float(self.bkgUnc[pr]))
      else: prstring.append('-')
    return self.FixListSpaces(prstring)

  def GetSystUncLine(self, systName):
    prstring = []
    for ipr in self.ProcessList():
      hnameUp = ipr + '_' + systName + 'Up'
      hnameDo = ipr + '_' + systName + 'Down'
      if hnameUp in self.f and hnameDo in self.f: prstring.append(1)
      else: prstring.append('-')
    return self.FixListSpaces(prstring)



if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Create a datacard from a rootfile')
  parser.add_argument('path',    help = 'Path to rootfile')
  parser.add_argument('--signal', '-s', default = 'tt' , help = 'Signal')
  parser.add_argument('--bkg',  '-b', default = [], help = 'List of bkgs')
  parser.add_argument('--bkgUnc',    '-k', default = None , help = 'Uncertainty on bkg normalization')
  parser.add_argument('--lumiUnc',    '-l', default = 1 , help = 'Uncertainty on lumi')
  parser.add_argument('--output',  '-o',  default = 'temp.dat', help = 'output file')
  parser.add_argument('--syst',  '-u', default = [] , help = 'List of systematics')
  args = parser.parse_args()

  fname = args.path
  signal = args.signal
  bkg = args.bkg
  norm = args.bkgUnc
  lumiUnc = float(args.lumiUnc)
  outname = args.output
  systList = args.syst
  d = Datacard(fname, signal, bkg, lumiUnc, norm, systList, outname=outname, nSpaces=12)
  d.Save()


