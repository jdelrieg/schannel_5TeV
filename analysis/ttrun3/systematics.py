from config import *
from cafea.plotter.xsec import *
from cafea.modules.fileReader import *
from PDFscaleUncertainties import *

names = {
  'lepSF_muon' : 'Muon efficiences',
  'lepSF_elec' : 'Electron efficiences',
  #'eleceff' : 'Electron efficiences',
  #'muoneff' : 'Muon efficiences',
  'trigSF' : 'Trigger efficiencies',
  'JES' : 'Jet energy scale',
  #'UE' : 'Underlying event',
  'hdamp' : 'ME/PS matching',#$h_\mathrm{damp}$',
  'mtop' : 'Top mass',
  'ISR' : 'Initial-state radiation',
  'FSR' : 'Final-state radiation',
  'DY' : 'Drell--Yan',
  'PU' : 'Pileup reweighting',
  'semilep' : r'$t\bar{t} \rightarrow 1 \ell$',
  #'JES_FlavorQCD': 'JES FlavorQCD',
  #'JES_SubTotalPileUp': 'JES Pileup',
  #'JES_SubTotalRelative':'JES Relative',
  #'JES_SubTotalAbsolute': 'JES Absolute',
  #'JES_TimePtEta': 'JES Time, pt and eta'
}

#path = 'histos/run3/5jun2022/'
### Fix categories
categories = {'channel':'em', 'level': 'g2jets', 'sign':'OS'}#, 'syst':'norm'}
categoriesPDF = {'channel':'em', 'level': 'g2jets'}#, 'syst':'norm'}
processDic = {
  'tt': 'TTTo2L2Nu',
  'tW': 'tbarW, tW',
  'semilep':'TTToSemiLeptonic',
  'WJets':'WJetsToLNu',
  'DY': 'DYJetsToLL_M50, DYJetsToLL_M10to50', 
  'Diboson' : 'WW, WZ, ZZ',#'WWTo2L2Nu, WZTo3LNu',#'WW, WZ, ZZTo2L2Nu',
  'data' : 'MuonEG,EGamma,DoubleMuon,SingleMuon,Muon'
}
bkglist    = ['tt', 'tW','semilep', 'WJets', 'DY', 'Diboson']
### Create plotter
p = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi, var='counts')

### Add hdamp and tune uncertainties
hdampup,hdampdo = GetModSystHistos(path, 'TTTo2L2Nu_hdamp', 'hdamp', var='counts')
#mtopup,mtopdo = GetModSystHistos(path, 'TTTo2L2Nu_mtop', 'mtop', var='counts')
#tuneup , tunedo = GetModSystHistos(path, 'TTTo2L2Nu_UE', 'UE', var='counts')
p.AddExtraBkgHist([hdampup, hdampdo], add=True)


### Create xsec object
experimental = ['lepSF_muon', 'lepSF_elec','PU','trigSF']#,'JER','JES_FlavorQCD','JES_SubTotalPileUp','JES_SubTotalRelative','JES_SubTotalAbsolute','JES_TimePtEta']
modeling = ['ISR', 'FSR','hdamp'] # ['UE', 'hdamp', 'ISR', 'FSR']
x = xsec('tt', 0.06, {'tW':0.15,'semilep':0.2,'WJets':0.3, 'DY':0.2, 'Diboson':0.3}, plotter=p, verbose=4, thxsec=921, experimental=experimental, modeling=modeling, categories=categories)
x.SetNames(names)
x.ComputeCrossSection()

pdf   = Get1bPDFUnc(  path, categories=categoriesPDF, sample='TTTo2L2Nu', doPrint=False)
scale = Get1binScaleUnc(path, categories=categoriesPDF, sample='TTTo2L2Nu', doPrint=False)
x.AddModUnc('PDF$+\\alpha_{S}$', pdf, isRelative=True)
x.AddModUnc('$\mu_R, \mu_F$ scales', scale, isRelative=True)
#x.AddModUnc('mtop', 0.005972, isRelative=True)

#jecs = GetJECSystHistos(path, 'variations/TTTo2L2Nu_withJEC', var='counts', categories=categories)

#x.AddExpUnc('Jet energy scale (external)', jecs, isRelative=True)
#x.AddExpUnc('Jet energy scale', 0.02, isRelative=True)
x.ComputeXsecUncertainties()

'''
path_gen = '/nfs/fanae/user/juanr/cafea/histos/TTTo2L2Nu_gen.pkl.gz'
hcounts = GetHisto(path_gen, "counts", categories={"sample":"TTTo2L2Nu", "channel":"all", "level":"all", "syst":"norm", "sign":"all"})
genEvents= list(hcounts.values().values())[0][0]
hfidu = GetHisto(path_gen, "counts", categories={"sample":"TTTo2L2Nu", "channel":"em", "level":"g2jets", "syst":"norm", "sign":"OS"})
fiduEvents = list(hfidu.values().values())[0][0]
BRem = 0.031938
BRSF = 0.015969
BRdilep = 0.108*0.108*9

x.SetGenEvents(genEvents)
x.SetFiduEvents(fiduEvents)
x.SetBR(BRem)
x.SetBRsample(BRdilep)

acc, accunc = x.GetAcceptance()
eff, effunc = x.GetEfficiency()
x.ComputeFiducial()
'''
x.GetYieldsTable()
x.GetUncTable(form='%1.2f')

