from config import *
#from PDFscaleUncertainties import *

from cafea.modules.CreateDatacardFromRootfile import Datacard
from cafea.plotter.plotter import GetHisto

if  not '/' in inputFile: outpath = './'
else: outpath = inputFile[:inputFile.rfind('/')]

from PDFscaleUncertainties import Get1binScaleUnc, Get1bPDFUnc
doPDFsplit = True
doScalesplit = False



channels = ['e', 'm']
levels = ['3j1b', '3j2b', '4j1b', '4j2b', 'g5j1b', 'g5j2b']
def GetChanLevFromName(fname):
  # medianDRjj_e_g5j1b.root
  inputs = fname[fname.rfind('/')+1:].replace('.root', '').split('_')
  chan = None; lev = None
  for i in inputs:
    if i in channels: chan = i
    if i in levels: lev = i
  if chan is None or lev is None:
    print("WARNING: could not get channel or level from file name: %s"%fname)
  return chan, lev

def GetModUnc(path, chan, lev):
  nbin = channels.index(chan) *len(levels) + levels.index(lev)
  if os.path.isfile(path + 'masterhistos/master.pkl.gz'):
    # print('Loading masterhistos from %s'%path)
    histo = GetHisto(path + 'masterhistos/master.pkl.gz', 'master').integrate('process', 'tt')
    nominal = histo.integrate('syst', 'norm').values(overflow='all')[()][nbin]
    # PDF and scale
    pdfUp = histo.integrate('syst', 'PDFUp').values(overflow='all')[()][nbin]
    pdfDown = histo.integrate('syst', 'PDFDown').values(overflow='all')[()][nbin]
    pdf = (abs(pdfUp-nominal) + abs(pdfDown-nominal))/2/nominal
    scaleUp = histo.integrate('syst', 'ScaleUp').values(overflow='all')[()][nbin]
    scaleDown = histo.integrate('syst', 'ScaleDown').values(overflow='all')[()][nbin]
    scales = (abs(scaleUp-nominal) + abs(scaleDown-nominal))/2/nominal
    # hdamp and UE
    tot = sum(histo.integrate('syst', 'norm').values(overflow='all')[()])
    hdampUp = sum(histo.integrate('syst', 'hdampUp').values(overflow='all')[()])
    hdampDown = sum(histo.integrate('syst', 'hdampDown').values(overflow='all')[()])
    hdamp = max(abs(hdampUp-tot), abs(hdampDown-tot))/tot
    UEUp = sum(histo.integrate('syst', 'UEUp').values(overflow='all')[()])
    UEDown = sum(histo.integrate('syst', 'UEDown').values(overflow='all')[()])
    UE = max(abs(UEUp-tot),abs(UEDown-tot))/tot
    return pdf, scales, hdamp, UE
  else:
    print("WARNING: please provide master histograms to take modeling uncertaintes... for now, returning hardcoded values")
    pdf = 0.007
    scales = 0.002
    hdamp = 0.007
    UE = 0.005
  return pdf, scales, hdamp, UE


def CreateDatacard(fname, outpath=outpath, oname=output):
  chan, lev = GetChanLevFromName(fname)
  if oname is None:
    oname = fname[fname.rfind('/')+1:] if '/' in fname else fname
    if oname.endswith('.root'): oname = oname[:-5]
    if '/' in oname: oname[oname.rfind('/')+1:]
  oname = 'dat_'+oname
  if not oname.endswith('.txt'): oname += '.txt'
  
  lumiUnc = 0.019
  bkg =  ['tW', 'tchan', 'WJets', 'QCD', 'DY']#'tchan',
  norm = [0.056, 0.02, 0.2, 0.3, 0.2]
  signal = 'tt'
  systList = ['muonSF', 'elecSF', 'btagSF','trigSF', 'FSR', 'ISR', 'prefire']# 'hdamp', 'UE', 'trigSF', 'Scales', 'PDF', 'Prefire']
  systList += ['JER','MC', 'AbsStat', 'AbsScale', 'AbsMPF', 'Frag', 'ECAL', 'HCAL', 'Flavor', 'RelStat', 'RelPt', 'RelBal', 'RelJER', 'L3Res','MET_UnclusteredEnergy']
  d = Datacard(fname, signal, bkg, lumiUnc, norm, systList, nSpaces=12, verbose=verbose)
  #d.AddExtraUnc('prefiring', 0.014, ['tt', 'tW', 'WJets', 'DY'])
  
  #pdf   = Get1bPDFUnc(  fname, categories=categoriesPDF, sample='TTTo2L2Nu', doPrint=False)
  #scale = Get1binScaleUnc(fname, categories=categoriesPDF, sample='TTTo2L2Nu', doPrint=False)
  
  pdf, scales, hdamp, UE = GetModUnc(path, chan, lev)
  if not doPDFsplit:
    d.AddExtraUnc('PDF', pdf, signal)
  else:
    ttSampleName = 'TTPS/' if os.path.isdir(path + 'TTPS/') else 'TTPS'
    pdfs = Get1bPDFUnc(path+ttSampleName, categories={'sample':processDic[signal], 'channel':chan, 'level':lev}, doPrint=False, returnAll=True)
    for i in range(len(pdfs)):
      d.AddExtraUnc('PDF%d'%(i+1), pdfs[i], signal)
  if not doScalesplit:
    d.AddExtraUnc('Scales', scales, signal)
  else:
    ttSampleName = 'TTPS/' if os.path.isdir(path + 'TTPS/') else 'TTPS'
    scales = Get1binScaleUnc(path+ttSampleName, categories={'sample':processDic[signal], 'channel':chan, 'level':lev}, doPrint=False, returnAll=True)
    for i in range(len(scales)):
      d.AddExtraUnc('Scales%d'%(i+1), scales[i], signal)
   
  #if chan == 'e': d.AddExtraUnc('trigElecSF', 0.02, signal) quitar
  #else          : d.AddExtraUnc('trigMuonSF', 0.01, signal)
  d.AddExtraUnc('hdamp', hdamp, signal)
  d.AddExtraUnc('UE', UE, signal)
  #d.AddExtraUnc('Prefire', 0.01, signal)
  d.SetOutPath(outpath)
  print('outpath = %s'%outpath)
  d.Save(oname)

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
