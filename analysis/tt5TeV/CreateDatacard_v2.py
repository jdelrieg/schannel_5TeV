#!/usr/bin/env python3

from config import *
#from PDFscaleUncertainties import *
import ROOT as r
from cafea.modules.CreateDatacardFromRootfile_v2 import Datacard
from cafea.plotter.plotter import GetHisto
from math import isnan

if  not '/' in inputFile: outpath = './'
else: outpath = inputFile[:inputFile.rfind('/')]

from PDFscaleUncertainties import Get1binScaleUnc, Get1bPDFUnc
doPDFsplit = True
doPDFsplit = True
doScalesplit = False

minBkgContrib = 1e-2
#minBkgContrib = 1e-1
minShapeVar   = 1e-5
minNormUnc    = 1e-3
removezeros   = True
reviewshapes  = True

#channels = ['e', 'm']
channels = ['e_plus','e_minus','m_plus','m_minus']
levels = ['2j1b','3j1b','3j2b']#,'2j0b']
systList = ['muonSF', 'elecSF', 'btagSFlight','btagSFbc','trigSF',  'prefire', 'JER', 'MC',
            'AbsStat', 'AbsScale', 'AbsMPF', 'Frag', 'ECAL', 'HCAL', 'Flavor', 'RelStat',
            'RelPt', 'RelBal', 'RelJER', 'L3Res','MET_UnclusteredEnergy','FSR', 'ISR','hdamp','UE','QCD_shape']#,'QCD_ll','QCD_lh','QCD_hl','QCD_hh'] #'btagSFlight','btagSFbc',
#systList = ['muonSF', 'elecSF', 'btagSF','trigSF', 'FSR', 'ISR', 'prefire', 'JER', 'Total','MET_UnclusteredEnergy']
smoothDict = {}
for iU in systList:
    smoothDict[iU] = {
        "e" : {'3j1b'  : {}, '3j2b'  : {}, '4j1b'  : {}, '4j2b'  : {}, 'g5j1b' : {}, 'g5j2b' : {},},
        "m" : {'3j1b'  : {}, '3j2b'  : {}, '4j1b'  : {}, '4j2b'  : {}, 'g5j1b' : {}, 'g5j2b' : {},}
    }

smoothDict = None # Sin suavizaciones

### Legend:
# Norm. smoothing:   ["norm", SYMM or NOT], e.g. ["norm", 1], ["norm", 0]
# Fitting smoothing: ["fit", SYMM OR NOT, ORDER], e.g. ["fit", 1, 2], ["fit", 0, 1]


#'muonSF'
#smoothDict["muonSF"]["e"] = {
#    '3j1b'  : ["norm", 1],
#    '3j2b'  : {},
#    '4j1b'  : {},
#    '4j2b'  : {},
#    'g5j1b' : {},
#    'g5j2b' : {},
#}
#smoothDict["muonSF"]["m"] = {
#    '3j1b'  : {},
#    '3j2b'  : {},
#    '4j1b'  : {},
#    '4j2b'  : {},
#    'g5j1b' : {},
#    'g5j2b' : {},
#}

#'elecSF'

#'btagSF'

#'trigSF'

#'FSR'

#'ISR'

#'prefire'

#'JER'

#'MC'

#'AbsStat'

#'AbsScale'

#'AbsMPF'

#'Frag'

#'ECAL'

#'HCAL'

#'Flavor'

#'RelStat'

#'RelPt'

#'RelBal'

#'RelJER'

#'L3Res'
#smoothDict["L3Res"]["e"] = {
#    '3j1b'  : {"WJets" : ["norm", 1]},
#    '3j2b'  : {},
#    '4j1b'  : {},
#    '4j2b'  : {},
#    'g5j1b' : {},
#    'g5j2b' : {},
#}
#smoothDict["L3Res"]["m"] = {
#    '3j1b'  : {},
#    '3j2b'  : {},
#    '4j1b'  : {},
#    '4j2b'  : {},
#    'g5j1b' : {},
#    'g5j2b' : {},
#}

#'MET_UnclusteredEnergy'


##############################################################################################
def GetChanLevFromName(fname):
  # medianDRjj_e_g5j1b.root
  inputs = fname[fname.rfind('/')+1:].replace('.root', '').split('_')
  chan = None; lev = None
  for i in inputs:
    if i in ['plus','minus']: i=inputs[-3]+'_'+i
    if i in channels: chan = i
    if i in levels: lev = i
  if chan is None or lev is None:
    print("WARNING: could not get channel or level from file name: %s"%fname)
  return chan, lev


def GetModUnc(path, chan, lev):  #This way calculated as flat (different for each level in pdfs and scales and same for all levels in ue hdamp. This comes from how are they calculated in MasterPlot)
  nbin = channels.index(chan) *len(levels) + levels.index(lev)
  if os.path.isfile(path + 'masterhistos/master.pkl.gz'):
    # print('Loading masterhistos from %s'%path)
    histo_top = GetHisto(path + 'masterhistos/master.pkl.gz', 'master').integrate('process', 'tchan')
    histo_tbar = GetHisto(path + 'masterhistos/master.pkl.gz', 'master').integrate('process', 'tbarchan')
    nominal_top = histo_top.integrate('syst', 'norm').values()[()][nbin]#overflow='all'						OJO, en todos los sitios que tengo el overflow comentado, estaba dentro de values(). Lo he modificado pq creo que eso estaba 
    nominal_tbar = histo_tbar.integrate('syst', 'norm').values()[()][nbin]#overflow='all'					moviendo un indice lo que teniamos en el master. De ese modo, el indice 'nbin' no se correspondia con el canal real. y de hecho petaba en el primero, pq el overflow añadia un 0 en el primer y ultimo indice
    # PDF and scale																					
    pdfUp_top = histo_top.integrate('syst', 'PDFUp').values(overflow='all')[()][nbin]
    pdfDown_top = histo_top.integrate('syst', 'PDFDown').values(overflow='all')[()][nbin]
    pdf_top = (abs(pdfUp_top-nominal_top) + abs(pdfDown_top-nominal_top))/2/nominal_top
    
    pdfUp_tbar = histo_tbar.integrate('syst', 'PDFUp').values(overflow='all')[()][nbin]
    pdfDown_tbar = histo_tbar.integrate('syst', 'PDFDown').values(overflow='all')[()][nbin]
    pdf_tbar = (abs(pdfUp_tbar-nominal_tbar) + abs(pdfDown_tbar-nominal_tbar))/2/nominal_tbar
    
    scaleUp_top = histo_top.integrate('syst', 'ScaleUp').values()[()][nbin] #overflow='all'
    scaleDown_top = histo_top.integrate('syst', 'ScaleDown').values()[()][nbin]#overflow='all'
    scales_top = (abs(scaleUp_top-nominal_top) + abs(scaleDown_top-nominal_top))/2/nominal_top

    scaleUp_tbar = histo_tbar.integrate('syst', 'ScaleUp').values()[()][nbin] #overflow='all'
    scaleDown_tbar = histo_tbar.integrate('syst', 'ScaleDown').values()[()][nbin]#overflow='all'
    scales_tbar = (abs(scaleUp_tbar-nominal_tbar) + abs(scaleDown_tbar-nominal_tbar))/2/nominal_tbar
    # hdamp and UE

#    tot = histo.integrate('syst', 'norm').values(overflow='all')[()][nbin+1]
#    hdampUp = histo.integrate('syst', 'hdampUp').values(overflow='all')[()][nbin+1]
#    hdampDown = histo.integrate('syst', 'hdampDown').values(overflow='all')[()][nbin+1]
#    hdamp = max(abs(hdampUp-tot), abs(hdampDown-tot))/tot
#    UEUp = histo.integrate('syst', 'UEUp').values(overflow='all')[()][nbin+1]
#    UEDown = histo.integrate('syst', 'UEDown').values(overflow='all')[()][nbin+1]
#    UE = max(abs(UEUp-tot),abs(UEDown-tot))/tot
     
    hdamp  =  0.0001
    UE     = 0.0001    
    return pdf_top, scales_top,pdf_tbar,scales_tbar#, hdamp, UE
  else:
    print("WARNING: please provide master histograms to take modeling uncertaintes... for now, returning hardcoded values")
    pdf    =  0.007
    scales =  0.002
    hdamp  =  0.007
    UE     = 0.005
  return pdf, scales, hdamp, UE


def reviewBkgContrib(iF, blist, nlist, ch, lv):
  oblist = []; onlist = []
  thef = r.TFile(iF, "READ")
  for iP in blist:
    if iP=='QCD':oblist.append(iP);onlist.append(nlist[blist.index(iP)]);continue
    if thef.Get(iP).Integral() > minBkgContrib:
      oblist.append(iP)
      onlist.append(nlist[blist.index(iP)])
      #print(iF,iP,thef.Get(iP).Integral())
    #elif iP=="QCD":			#esto lo añadi con el truquillo pa quitarm qcd
     # continue				#same
    else:
      print("\t- WARNING: {}-{} background process {} dropped because of low contribution:{}.".format(ch, lv, iP,thef.Get(iP).Integral()))#print(f"\t- WARNING: {ch}-{lv} background process {iP} dropped because of low contribution.")
  thef.Close(); del thef
  return oblist, onlist

def reviewBkgContrib_2(iF, blist, nlist, ch, lv):    #NOTE: this check works assuming symmetric uncertainties (they more or less are)
													 # and it is not perfect since we are not including all uncs, only the 'shape-type' ones. TO compensate this effect, 
													 # we lower the threshold (ratioup>1->ratioup>0.5)
													 # So this is not a perfect check but at least gives an idea and it is good to have a look at it when producing cards :)
  oblist = []; onlist = []
  uncs=['btagSFlight','btagSFbc','elecSF', 'muonSF', 'prefire','JER','MC', 'AbsStat', 'AbsScale', 'AbsMPF', 'Frag', 'ECAL', 'HCAL', 'Flavor', 'RelStat', 'RelPt', 'RelBal', 'RelJER', 'L3Res','MET_UnclusteredEnergy','trigSF']
  thef = r.TFile(iF, "READ") #'btagSFlight','btagSFbc',
  for iP in blist:
    if iP=='QCD':oblist.append(iP);onlist.append(nlist[blist.index(iP)]);continue
    nominal=thef.Get(iP).Integral()
    countup=0
    ratioup=0
    for uncer in uncs:
      countup=countup+(thef.Get(iP+"_" + uncer + "Up").Integral()-nominal)
    if nominal>0: ratioup=(countup-nominal)/nominal
    if ratioup >0.5:  
      print("\t- WARNING: {}-{} background process {} dropped because of high ratio:{}.".format(ch, lv, iP,ratioup))
    else:
      oblist.append(iP)
      onlist.append(nlist[blist.index(iP)])				
    
  thef.Close(); del thef
  return oblist, onlist

def reviewQCDshapes(iF, blist, nlist, ch, lv):
  oblist = []; onlist = []
  thef = r.TFile(iF, "READ")
  uncs=['QCD_ll','QCD_lh','QCD_hl','QCD_hh']
  for iP in blist:
    nominal=thef.Get(iP).Integral()
    if iP != "QCD":
      oblist.append(iP)
      onlist.append(nlist[blist.index(iP)])		
    else:
      for uncer in uncs:
        up=0
        up=thef.Get(iP+"_" + uncer + "Up").Integral()
        if abs(nominal-up)/nominal>1.1:
          print(" YOU MUST REMOVE IN: {}-{} QCD because of the uncertainty:{}. IT IS {}.".format(ch, lv,uncer,abs(nominal-up)/nominal))
			  
      oblist.append(iP)
      onlist.append(nlist[blist.index(iP)])
      #print(iF,iP,thef.Get(iP).Integral())
    #elif iP=="QCD":			#esto lo añadi con el truquillo pa quitarm qcd
     # continue				#same
  thef.Close(); del thef
  return oblist, onlist

def CreateDatacard(fname, outpath=outpath, oname=output):
#  print("Iniciando lo de las cards")
  chan, lev = GetChanLevFromName(fname)
  if oname is None:
    oname = fname[fname.rfind('/')+1:] if '/' in fname else fname
    if oname.endswith('.root'): oname = oname[:-5]
    if '/' in oname: oname[oname.rfind('/')+1:]
  oname = 'dat_' + oname
  if not oname.endswith('.txt'): oname += '.txt'
  
  lumiUnc  = 0.0#0.019        #If wanting to add it, go to CreateDatacradFromRootfile_v2 line 278 (AddLumiUncTxt) and uncomment it
  bkg      = ['tW', "tt", 'WJetsH','WJetsL', 'QCD', 'DY']
  norm     = [0.056, 0.05,   0.2, 0.2,  0.3,   0.3] #QCD era 0.3        #ESTO LO CAMBIE YO
  signal   = ['tchan','tbarchan']
  
  bkg, norm = reviewBkgContrib(fname, bkg, norm, chan, lev)
  bkg, norm = reviewBkgContrib_2(fname, bkg, norm, chan, lev)
  #bkg, norm = reviewQCDshapes(fname, bkg, norm, chan, lev)
 
 
  d = Datacard(fname, signal, bkg, lumiUnc, norm, systList, nSpaces = 12, verbose = verbose, minShapeVar = minShapeVar, 
               minNormUnc = minNormUnc, rmNegBins = removezeros, reviewshapes = reviewshapes, smoothingDict = smoothDict,
               ch = chan, lv = lev)
 
  pdf_top, scales_top,pdf_tbar,scales_tbar = GetModUnc(path, chan, lev) #, hdamp, UE

  if not doPDFsplit:
    if isnan(pdf_top):
      print("\t- WARNING: while assuming PDF as having only norm. effect in {}-{} region it is found to be NaN! THIS SHOULD NOT HAPPEN AND WILL RESULT IN ISSUES IN COMBINE!".format(chan, lev))#print(f"\t- WARNING: while assuming PDF as having only norm. effect in {chan}-{lev} region it is found to be NaN! THIS SHOULD NOT HAPPEN AND WILL RESULT IN ISSUES IN COMBINE!")
    d.AddExtraUnc('PDF_top', pdf_top, signal[0])
    d.AddExtraUnc('PDF_tbar', pdf_tbar, signal[1])
  else:
    ttSampleName = 'TTPS/' if os.path.isdir(path + 'TTPS/') else 'TTPS'
    #pdfs = Get1bPDFUnc(path + ttSampleName, categories={'sample':[processDic[s] for s in signal], 'channel':chan, 'level':lev}, doPrint=False, returnAll=True)
    pdfs = Get1bPDFUnc(path + ttSampleName, categories={'sample':[processDic['tchan']] , 'channel':chan, 'level':lev}, doPrint=False, returnAll=True)
    for i in range(len(pdfs)):
      d.AddExtraUnc('PDF_top%d'%(i+1), pdfs[i], signal[0])  
    pdfs = Get1bPDFUnc(path + ttSampleName, categories={'sample':[processDic['tbarchan']] , 'channel':chan, 'level':lev}, doPrint=False, returnAll=True)
    for i in range(len(pdfs)):
      d.AddExtraUnc('PDF_tbar%d'%(i+1), pdfs[i], signal[1])
      
  if not doScalesplit:
    if isnan(scales_top):
      print("\t- WARNING: while assuming renor. and fact. scales as having only norm. effect in {}-{} region it is found to be NaN! THIS SHOULD NOT HAPPEN AND WILL RESULT IN ISSUES IN COMBINE!".format(chan,lev))
      if chan == "e" and lev == "3j1b":
        print("- FIXME: setting for e-3j1b scales value to provided one.")
        scales = 0.0010787972982649
    d.AddExtraUnc('Scales_top', scales_top, signal[0])
    d.AddExtraUnc('Scales_tbar', scales_tbar, signal[1])
  else:
    ttSampleName = 'TTPS/' if os.path.isdir(path + 'TTPS/') else 'TTPS'
    #scales = Get1binScaleUnc(path+ttSampleName, categories={'sample':processDic[signal], 'channel':chan, 'level':lev}, doPrint=False, returnAll=True)
    scales = Get1binScaleUnc(path+ttSampleName, categories={'sample':processDic['tchan'], 'channel':chan, 'level':lev}, doPrint=False, returnAll=True)
    for i in range(len(scales)):
      d.AddExtraUnc('Scales_top%d'%(i+1), scales[i], signal[0])
    scales = Get1binScaleUnc(path+ttSampleName, categories={'sample':processDic['tbarchan'], 'channel':chan, 'level':lev}, doPrint=False, returnAll=True)
    for i in range(len(scales)):
      d.AddExtraUnc('Scales_tbar%d'%(i+1), scales[i], signal[1])

#  if isnan(hdamp):
#    print("\t- WARNING: while assuming hdamp unc. for ttbar as having only norm. effect in {}-{} region it is found to be NaN! THIS SHOULD NOT HAPPEN AND WILL RESULT IN ISSUES IN COMBINE!".format(chan,lev))
#  d.AddExtraUnc('hdamp', hdamp, signal)
#  if isnan(UE):
#    print("\t- WARNING: while assuming UE unc. for ttbar as having only norm. effect in {}-{} region it is found to be NaN! THIS SHOULD NOT HAPPEN AND WILL RESULT IN ISSUES IN COMBINE!".format(chan,lev))
#  d.AddExtraUnc('UE', UE, signal)
  
  d.SetOutPath(outpath)
  d.Save(oname)
  return


if inputFile == '': 
  print('Please provide a root file to create datacards from using --inputFile /path/to/inputs/');
  exit(1)


if os.path.isdir(inputFile):
  for d in os.listdir(inputFile):
    if not d.endswith('.root'): continue
    fname = os.path.join(inputFile, d)
    CreateDatacard(fname)
else:
  CreateDatacard(inputFile)
