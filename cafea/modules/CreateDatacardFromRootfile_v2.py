# Create datacard from simple "combine" rootfile, using uproot
# Usage: 
#   python cafea/modules/CreateDatacardFromRootfile.py NBtags.root -s tt -b "tW, WJets, QCD, DY" -l 0.015 --bkgUnc "0.2, 0.2, 0.2, 0.2" -u "lepSF, trigSF, btagSF, FSR, ISR" -o test.dat
import uproot, os
import numpy as np
import awkward as ak
from coffea import hist, processor
from cafea.modules.FixNegValuesRoot_v2 import FixNegValuesRoot
from cafea.modules.applySmoothing      import applySmoothing
import ROOT as r

class Datacard:
  def __init__(self, fname, signal=[], bkgList=[], lumiUnc=0.0, bkgUnc=[], systList=[], nSpaces=10, outname=None, verbose = True,
               rmNegBins = True, minShapeVar = 1e-5, minNormUnc = 1e-3, tolSymmNorm = 1e-3, defaultZeroVal = 1e-5, tolstatfl = 0.5,
               tolcheck = 0.25, ch = "", lv = "", smoothingDict = None, reviewshapes = True):
    if rmNegBins:
      FixNegValuesRoot(fname, defaultZeroVal, verbose = verbose)
    if smoothingDict:
      applySmoothing(fname, smoothingDict, ch, lv, defaultZeroVal, verbose = verbose)
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
    self.doReview    = reviewshapes
    self.minShapeVar = minShapeVar
    self.minNormUnc  = minNormUnc
    self.tolstatfl   = tolstatfl
    self.tolcheck   = tolcheck
    #self.nBinsReal = self.f.get(signal if signal != "" else bkgList[0]).to_pyroot().GetNbinsX()
    self.nBinsReal = self.f.get(signal[0]).to_pyroot().GetNbinsX()#
    self.ch = ch
    self.lv = lv
    return

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
    if isinstance(signal, str) and ',' in signal: signal = signal.replace(' ', '').split(',')
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
    return len(self.bkg)+2#1

  def ProcessList_or(self):
    return [self.signal] + self.bkg
 
  def ProcessList(self):
    original =[self.signal] + self.bkg
    flat_list=[]
    for item in original:
      if isinstance(item,list):
        flat_list.extend(item)
      else:
        flat_list.append(item)
    return flat_list

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
    self.AddLine('jmax %i processes minus 2'%(len(self.bkg)+1))
    self.AddLine('kmax * number of nuisance parameters')
    self.AddSep()
    self.AddLine('shapes * %s %s $PROCESS $PROCESS_$SYSTEMATIC'%(self.chName, self.fname.split("/")[-1]))
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
    outline = ["Lumi", "lnN"] + [1 + self.lumiUnc]
    for iP in self.bkg:
        outline.append((1 + self.lumiUnc) if iP != "QCD" else "-")
    datacard = self.AddLine(self.FixListSpaces(outline))
    return


  def AddBkgNormTxt(self):
    for pr in self.bkg:
      self.AddLine(self.FixListSpaces([pr, 'lnN']) + self.GetNormUncProcess(pr))


  def isShapeVar(self, hu, hd, hn):
    ratioU = None
    ratioD = None
    for iB in range(1, self.nBinsReal + 1):
      if ratioU:
        if abs(hu.GetBinContent(iB)/hn.GetBinContent(iB) - ratioU) > self.minShapeVar:
          return True
      if ratioD:
        if abs(hd.GetBinContent(iB)/hn.GetBinContent(iB) - ratioD) > self.minShapeVar:
          return True
      ratioU = hu.GetBinContent(iB)/hn.GetBinContent(iB)
      ratioD = hd.GetBinContent(iB)/hn.GetBinContent(iB)
    return False


  def doUncCheck(self, hu, hd, hn, sy, pr):
    blonesdd = []
    blstatfl = []
    for iB in range(1, self.nBinsReal + 1):
      if (hu.GetBinContent(iB)/hn.GetBinContent(iB) - 1) * (hd.GetBinContent(iB)/hn.GetBinContent(iB) - 1) > 0:
        blonesdd.append(iB)
#        print(hu.GetBinContent(iB), hu.GetBinError(iB), hd.GetBinContent(iB), hd.GetBinError(iB))
      if (hu.GetBinError(iB) > self.tolstatfl*hu.GetBinContent(iB)):
        blstatfl.append(iB)
      if (hd.GetBinError(iB) > self.tolstatfl*hd.GetBinContent(iB)):
        blstatfl.append(iB)
    
    if float(len(blonesdd)) > self.tolcheck * self.nBinsReal:
      print(f"\t- WARNING: {self.ch}-{self.lv} unc. {sy} in proc. {pr} has significant one-sided variations.")
    if float(len(blstatfl)) > self.tolcheck * self.nBinsReal * 2:
      print(f"\t- WARNING: {self.ch}-{self.lv} unc. {sy} in proc. {pr} has significant statistical uncertainties in their variations.")
#    else:
#      print(len(blstatfl))
    
    return
    

  def AddSystUncTxt(self):
    ### TODO: implement automatic change to normalisation in case for all processes the actual impact is only on that.
    ### (it is almost implemented, but not enabled)
    for syst in self.syst:
      prstring = []
      isShape  = False
      for ipr in self.ProcessList():
        if ipr == "QCD" and syst not in ['QCD_ll','QCD_lh','QCD_hl','QCD_hh','QCD_shape']:
          prstring.append("-")
          continue       
        hnameUp = ipr + '_' + syst + 'Up'
        hnameDo = ipr + '_' + syst + 'Down'
        if hnameUp in self.f and hnameDo in self.f:
          if self.nBinsReal < 2:
            print("WARNING: this should be a norm. unc., but for consistency across channels will remain a shape unc..")
            prstring.append(1)
          elif self.doReview:
            if not isShape:# Check whether there is variation in shape in this process
              hup = self.f.get(hnameUp).to_pyroot()
              hdn = self.f.get(hnameDo).to_pyroot()
              h0  = self.f.get(ipr).to_pyroot()
              
              if self.isShapeVar(hup, hdn, h0):
                isShape = True
                
            if isShape:
              prstring.append(1)
              self.doUncCheck(hup, hdn, h0, syst, ipr)
            else:
              uprat = hup.Integral()/h0.Integral()
              dnrat = hdn.Integral()/h0.Integral()
              
              if max([abs(uprat - 1), abs(dnrat - 1)]) < self.minNormUnc:
                print(f"\t- WARNING: {self.ch}-{self.lv} unc. {syst} in proc. {ipr} does not have a genuine shape effect, and its normalisation effect is very small. We will ignore this uncertainty (in this process, channel, and level).")
                prstring.append("-")
              else:
                print(f"\t- WARNING: {self.ch}-{self.lv} unc. {syst} in proc. {ipr} does not have a genuine shape effect. WE WILL CONSIDER IT SHAPE ANYWAY.")
                prstring.append(1)
                
                self.doUncCheck(hup, hdn, h0, syst, ipr)
                
                ###############################################################################
#                if abs(uprat - 1/dnrat) < tolSymmNorm:
#                  prstring.append(str(uprat if uprat > 1 else dnrat))
#                else:
#                  prstring.append(f"{uprat}/{dnrat}")
          else:
            prstring.append(1)
        else:
          prstring.append('-')
      if any([el != "-" for el in prstring]):
        if all([el != "1" for el in prstring]):
#          self.AddLine(self.FixListSpaces([syst, 'lnN']) + self.FixListSpaces(prstring)) ################ FORCED
          self.AddLine(self.FixListSpaces([syst, 'shape']) + self.FixListSpaces(prstring))
        else:
          self.AddLine(self.FixListSpaces([syst, 'shape']) + self.FixListSpaces(prstring))
    return

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
#    self.AddLine('%s autoMCStats 0 0 1'%self.chName)
    self.AddLine('%s autoMCStats 10 0 1'%self.chName)

  def GetDatacard(self):
    self.datacard = ''
    self.AddTxtHeader()
    self.AddTxtObs()
    self.AddTxtRates()
    #self.AddLumiUncTxt()
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
    print('processList',ak.flatten(self.ProcessList()))
    return ['%1.5f'%(self.Yield(pr)) for pr in self.ProcessList()]

  def GetNormUncProcess(self, pr):
    prstring = []
    for ipr in self.ProcessList():
      if ipr == pr: prstring.append(1+float(self.bkgUnc[pr]))
      else: prstring.append('-')
    return self.FixListSpaces(prstring)



if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Create a datacard from a rootfile')
  parser.add_argument('path',                                  help = 'Path to rootfile')
  parser.add_argument('--signal',  '-s', default = 'tt',       help = 'Signal')
  parser.add_argument('--bkg',     '-b', default = [],         help = 'List of bkgs')
  parser.add_argument('--bkgUnc',  '-k', default = None ,      help = 'Uncertainty on bkg normalization')
  parser.add_argument('--lumiUnc', '-l', default = 1 ,         help = 'Uncertainty on lumi')
  parser.add_argument('--output',  '-o', default = 'temp.dat', help = 'output file')
  parser.add_argument('--syst',    '-u', default = [] ,        help = 'List of systematics')
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


