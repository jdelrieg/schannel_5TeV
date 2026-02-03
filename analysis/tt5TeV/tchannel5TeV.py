#!/usr/bin/env python
import lz4.frame as lz4f
import cloudpickle
import json, pickle
import pprint
import numpy as np
import awkward as ak
import coffea
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea import hist, processor
from coffea.util import load, save
from optparse import OptionParser
from coffea.analysis_tools import PackedSelection
from coffea.lumi_tools import LumiMask

from cafea.analysis.objects import *
from cafea.analysis.corrections import GetBTagSF, GetBtagEff, AttachMuonSF, AttachElectronSF, GetPUSF, GetTriggerSF5TeV, GetElecScale5TeV, jet_factory, jet_factory_data, met_factory, GetBtagSF5TeV
from cafea.analysis.selection import *
from cafea.modules.paths import cafea_path


from analysis.tt5TeV.w_reco_for import reco


import warnings
# ignore userwarnings and runtime warnings
warnings.filterwarnings("ignore", category=UserWarning) # evaluating the MVA because the dataset has not names for the variables...
warnings.filterwarnings("ignore", category=RuntimeWarning) # variables in the nanoAOD 

splitJES = True
fillAll = True
fillMVA = False
doSyst = False
fillAcc = False
'''
def AttachTrigSF(e0, m0, events):
  TrigSFe, TrigSFedo, TrigSFeup = GetTriggerSF5TeV(e0.pt, np.abs(e0.eta), 'e')
  TrigSFm, TrigSFmdo, TrigSFmup = GetTriggerSF5TeV(m0.pt, np.abs(m0.eta), 'm')
  TrigSFe   = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFe, 1.)), nan=1)
  TrigSFedo = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFedo, 1.)), nan=1)
  TrigSFeup = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFeup, 1.)), nan=1)
  TrigSFm   = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFm, 1.)), nan=1)
  TrigSFmdo = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFmdo, 1.)), nan=1)
  TrigSFmup = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFmup, 1.)), nan=1)
  events['sf_trig']    = TrigSFe*TrigSFm
  events['sf_trig_hi'] = TrigSFeup*TrigSFmup
  events['sf_trig_lo'] = TrigSFedo*TrigSFmdo
'''
def AttachTrigSF(e0, m0, events):
  TrigSFe, TrigSFedo, TrigSFeup = GetTriggerSF5TeV(np.abs(e0.eta), 'e')
  TrigSFm, TrigSFmdo, TrigSFmup = GetTriggerSF5TeV(np.abs(m0.eta), 'm')
  TrigSFe   = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFe, 1.)), nan=1)
  TrigSFedo = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFedo, 1.)), nan=1)
  TrigSFeup = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFeup, 1.)), nan=1)
  TrigSFm   = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFm, 1.)), nan=1)
  TrigSFmdo = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFmdo, 1.)), nan=1)
  TrigSFmup = np.nan_to_num(ak.flatten(ak.fill_none(TrigSFmup, 1.)), nan=1)
  events['sf_trig']    = TrigSFe*TrigSFm
  events['sf_trig_hi'] = TrigSFeup*TrigSFmup
  events['sf_trig_lo'] = TrigSFedo*TrigSFmdo

def GetNBtagNJets(njets, nbtags):
  '''
  njets and nbtags are numpy arrays. Let's use these mapping depending on (njets, nbtags)
  njets goes up to 6, nbtags up to njets or 2 as maximum
  (0, 0) --> 0
  (1, 0) --> 1, (1, 1) --> 2
  (2, 0) --> 3, (2, 1) --> 4, (2, 2) --> 5
  (3, 0) --> 6, (3, 1) --> 7, (3, >=2) --> 8,
  (4, 0) --> 9, (4, 1) --> 10, (4, >=2) --> 11,
  (5, 0) --> 12, (5, 1) --> 13, (5, >=2) --> 14,
  (>=6, >=0) --> 15
  '''
  nbtags = np.minimum(nbtags, 2)
  vals = np.copy(njets)
  vals = np.where(vals>=2, 3*(njets-1) + nbtags, vals)
  vals = np.where(vals==1, njets + nbtags, vals)
  return vals






def GetJESMETvar(corrected_jets, met, jetptcut, btagwp, isData, var=''):
  if var.lower().endswith('up'):
    dir = 'up'
    var = var[:-2]
  elif var.lower().endswith('down'):
    dir = 'down'
    var = var[:-4]
  elif var.lower().endswith('do'):
    dir = 'down'
    var = var[:-2]
  cleanedJets = ak.copy(corrected_jets)
  jetptname = "pt_nom" if (hasattr(cleanedJets, "pt_nom") and not isData) else "pt"
  metpt = met.pt
  if var != '':
    if 'MET' in var:
      if dir.lower() == 'up':
        metpt       = getattr(met, var).up.pt
      else:
        metpt       = getattr(met, var).down.pt      
    elif var == 'JER':
      if dir.lower() == 'up':
        cleanedJets = getattr(corrected_jets, var).up
        metpt       = getattr(met, var).up.pt
      else:
        cleanedJets = getattr(corrected_jets, var).down
        metpt       = getattr(met, var).down.pt
    else:
      if dir.lower() == 'up':
        cleanedJets = getattr(corrected_jets, 'JES_'+var).up
        metpt       = getattr(met, 'JES_'+var).up.pt
      else:
        cleanedJets = getattr(corrected_jets, 'JES_'+var).down
        metpt       = getattr(met, 'JES_'+var).down.pt
  cleanedJets["isGood"] = isTightJet(getattr(cleanedJets, jetptname),cleanedJets.eta, cleanedJets.jetId, jetPtCut=jetptcut)#,jetEtaCut=jetetacut)    #in isTightJet we define the etacut
  goodJets = cleanedJets[cleanedJets.isGood]
  goodJets["isBtag"] = (goodJets.btagDeepB > btagwp)
  btagSF = np.ones(len(goodJets), dtype=np.float64)
  if not isData: # btag SF
    btagSF = GetBtagSF5TeV(goodJets.pt, goodJets.eta, goodJets.hadronFlavour, goodJets.isBtag, False)
  return goodJets, metpt, btagSF

def GetNjetNbtagsNujets(jets):
  njets = ak.num(jets)
  nbtags = ak.num(jets[(jets.isBtag== 1)])
  nujets = ak.num(jets[(jets.isBtag== 0)])
  return njets, nbtags, nujets


def GetCutJets(cuts, syst, metpt, njets, nbtags, nujets=None):
  jetcat = np.ones_like(metpt, dtype=bool)
  metcat = np.ones_like(metpt, dtype=bool)

  if 'metg40' in cuts:
    metcat = (metpt > 40)
  if 'metl40' in cuts:
    metcat = (metpt < 40)
  elif 'metg45' in cuts:
    metcat = (metpt > 45)
  elif 'metl45' in cuts:
    metcat = (metpt < 45)
  elif 'metg15' in cuts:
    metcat = (metpt > 15)
  elif 'metl15' in cuts:
    metcat = (metpt < 15)
  elif 'metl35' in cuts:
    metcat = (metpt < 35)
  elif 'metg35' in cuts:
    metcat = (metpt > 35)
  elif 'metg30' in cuts:
    metcat = (metpt > 30)
  elif 'metl30' in cuts:
    metcat = (metpt < 30)

  elif 'metl25' in cuts:
    metcat = (metpt < 25)
  elif 'metg25' in cuts:
    metcat = (metpt > 25)
  elif 'metg20' in cuts:
    metcat = (metpt > 20)
  elif 'metl20' in cuts:
    metcat = (metpt < 20)

#QCD ex
  elif 'metl29' in cuts:
    metcat = (metpt < 29)
  elif 'metl28' in cuts:
    metcat = (metpt < 28)
  elif 'metl27' in cuts:
    metcat = (metpt < 27)
  elif 'metl26' in cuts:
    metcat = (metpt < 26)
  elif 'metl24' in cuts:
    metcat = (metpt < 24)
  elif 'metl23' in cuts:
    metcat = (metpt < 23)
  elif 'metl22' in cuts:
    metcat = (metpt < 22)
  elif 'metl21' in cuts:
    metcat = (metpt < 21)  
#QCD ex    

  if 'g2jets' in cuts:
    jetcat = (njets >= 2) & (nbtags>=0)
  elif 'g3jets' in cuts:
    jetcat = (njets >= 3)
  elif 'g4jets' in cuts:
    jetcat = (njets >= 4)
  elif '0b' in cuts:  
    jetcat = (nbtags == 0)
  elif '1b' in cuts:
    jetcat = (nbtags == 1)
  elif '2b' in cuts:
    jetcat = (nbtags == 2)
  elif '2j0b' in cuts:
    jetcat = (njets == 2) & (nbtags == 0)
  elif '2j1b' in cuts:
    jetcat = (njets == 2) & (nbtags == 1)
  elif '3j1b' in cuts:
    jetcat = (njets == 3) & (nbtags == 1)
  elif '3j2b' in cuts:
    jetcat = (njets == 3) & (nbtags == 2)  ##OJO cambio con tt, aqui es exactamente, no mas de dos
  elif '4j1b' in cuts:
    jetcat = (njets == 4) & (nbtags == 1)
  elif '4j2b' in cuts:
    jetcat = (njets == 4) & (nbtags >= 2)
  elif 'g5j1b' in cuts:
    jetcat = (njets >= 5) & (nbtags == 1)
  elif 'g5j2b' in cuts:
    jetcat = (njets >= 5) & (nbtags >= 2)
  elif 'g2ujets' in cuts:
    jetcat = (nujets >= 2)
  elif '3j0b' in cuts:
    jetcat = (njets == 3) & (nbtags == 0)
  elif '4j0b' in cuts:
    jetcat = (njets == 4) & (nbtags == 0)
  elif 'g5j0b' in cuts:
    jetcat = (njets >= 5) & (nbtags == 0)

  return (jetcat) & (metcat)


class EventCollection:
    def __init__(self, lorentz_vectors):
        self.lorentz_vectors = lorentz_vectors

    def M(self):
        return np.array([lv.M() if lv is not None else 0 for lv in self.lorentz_vectors])

    def Px(self):
        return np.array([lv.Px() if lv is not None else 0 for lv in self.lorentz_vectors])

    def Py(self):
        return np.array([lv.Py() if lv is not None else 0 for lv in self.lorentz_vectors])

    def Pz(self):
        return np.array([lv.Pz() if lv is not None else 0 for lv in self.lorentz_vectors])

    def E(self):
        return np.array([lv.E() if lv is not None else 0 for lv in self.lorentz_vectors])
        
    def eta(self):
        return np.array([lv.Eta() if lv is not None else 0 for lv in self.lorentz_vectors])        
    
    def phi(self):
        return np.array([lv.Phi() if lv is not None else 0 for lv in self.lorentz_vectors])   



class AnalysisProcessor(processor.ProcessorABC):
  
    def __init__(self, samples, model):

        self._samples = samples
        self.model = model
        self.model_2j1b = self.model; self.model_3j2b = None
        if isinstance(model, list) and len(model) >= 2:
          
          
          self.model_pruned = model[0]
          self.model_relaxed_b10 = model[1]
          

        # Create the histograms
        # 'name' : hist.Hist("Ytitle", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat("syst", "syst"), hist.Bin("name", "X axis (GeV)", 20, 0, 100)),
        self._accumulator = processor.dict_accumulator({
        'dummy'      : hist.Hist("Dummy", hist.Cat("sample", "sample"), hist.Bin("dummy", "Number of events", 1, 0, 1)),
        'counts'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'l0pt'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("lep0pt",  "Leading lepton $p_{T}$ (GeV)", 10, 20, 120)),
        'PDF'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("PDF",     "Counts", 103, 0, 103)),
        'Scales'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("Scales",  "Counts", 9, 0, 9)),
        'l0eta'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("lep0eta", "Leading lepton $\eta$ ", 8, -2.5, 2.50)),
        'ept'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("ept",  "Electron $p_{T}$ (GeV)", 8, 20, 120)),
        'eeta'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("eeta", "Electron $\eta$ ", 8, -2.5, 2.50)),
        'mpt'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("mpt",  "Muon $p_{T}$ (GeV)", 8, 20, 120)),
        'meta'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("meta", "Muon $\eta$ ", 8, -2.5, 2.50)),
        'j0pt'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("j0pt",  "Leading jet $p_{T}$ (GeV)", 10, 0, 300)),
        'j0eta'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("j0eta", "Leading jet $\eta$ ", 8, -2.5, 2.50)),
        'invmass'    : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("invmass", "$m_{\ell\ell}$ (GeV) ", 20, 0, 300)),
        'invmass2'   : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("invmass2", "$m_{\ell\ell}$ (GeV) ", 30, 70, 110)),
        'invmass_bb' : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("invmass", "$m_{\ell\ell}$ (GeV) ", 30, 70, 110)),
        'invmass_be' : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("invmass", "$m_{\ell\ell}$ (GeV) ", 30, 70, 110)),
        'invmass_ee' : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("invmass", "$m_{\ell\ell}$ (GeV) ", 30, 70, 110)),
        'njets'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("njets",   "Jet multiplicity", 10, 0, 10)),
        'nbtags'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("nbtags",  "b-tag multiplicity", 4, 0, 4)),
        'met'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("met",     "MET (GeV)", 10, 0, 200)),
        'metnocut'   : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("metnocut",     "MET (GeV)", 10, 0, 200)),
        'ht'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("ht",      "H$_{T}$ (GeV)", 20, 50, 400)),
        'mt'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("mt",      "m$_{T}$ (GeV)", 10, 0, 150)),
        'mlb'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("mlb",     "m(l,b) (GeV)", 12, 0, 400)),
        'minDRjj'    : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("minDRjj", "min$\Delta$R(jj) ", 10, 0, 5)),
        'medianDRjj' : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("medianDRjj", "median$\Delta$R(jj) ", 15, 0, 4.5)),
        'mjj'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("mjj",     "m(jj)( (GeV)", 12, 0, 240)),
        'ptjj'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("ptjj",    "p$_{T}$(jj) (GeV)", 15, 0, 300)),
        'njetsnbtags': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("njetsnbtags","(j,b)", 16, 0, 16)),
        'njetsnbtags12': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("njetsnbtags12","(j,b)", 12, 4, 16)),
  
        'u0pt'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("u0pt",  "Leading u-jet $p_{T}$ (GeV)", 10, 0, 300)),
        'b0pt'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("b0pt",  "Leading b-jet $p_{T}$ (GeV)", 10, 0, 300)),         
        'u0eta'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("u0eta", "Leading u-jet $\eta$ ", 8, -2.5, 2.50)),
        'absu0eta'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("absu0eta", "Leading u-jet $|\eta|$ ", 10, 0, 5)),        
        'minDRuu'    : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("minDRuu", "min$\Delta$R(uu) ", 10, 0, 3)),
        'medianDRuu' : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("medianDRuu", "median$\Delta$R(uu) ", 15, 0, 4.5)),
        'muu'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("muu",     "m(uu)( (GeV)", 10, 0, 200)),
        'ptuu'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("ptuu",    "p$_{T}$(uu) (GeV)", 15, 0, 300)),

         # ptSumVecAll, ptSumVeclb, dRlb = GetJetLepVar(goodJets, leps)
        'st'           : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("st",      "S$_{T}$ (GeV)", 20, 0, 800)),
        'ptlb'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("ptlb",    "p$_{T}$($\ell$b) (GeV)", 10, 0, 300)),
        'sumallpt'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("sumallpt", "$\sum_\mathrm{j,\ell}\,\mathrm{p}_{T}$ (GeV)", 10, 0, 300)),
        'dRlb'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("dRlb", "min$\Delta$R($\ell$b) ", 10, 0, 5)),
        'MVAscore_pruned'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("MVAscore_pruned", "MVA score", 20,0,1 )),
        'MVAscore_b6'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("MVAscore_b6", "MVA score", 20,0,1 )),
        'MVAscore_b10'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("MVAscore_b10", "MVA score", 20,0,1 )),
        
        'MVAscore_tX'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("MVAscore_tX", "MVA score", 20,0,1 )),
        'MVAscore_relaxed_b10'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("MVAscore_relaxed_b10", "MVA score", 20,0,1 )),
    

        'MVAtth'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("MVAtth", "MVA tth", 20, -1, 1)),
        'MVAwp'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("MVAwp", "MVA wp", 2, 0, 2)),


        'counts_metg25': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_metl25': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_metg15': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_metl15': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_metg30': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_metl30': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_metg20': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_metl20': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
     
        
                

        'pttrig' : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Bin("pttrig",  "Lepton $p_\mathrm{T}$", 20, 0, 100)),
        'etatrig': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Bin("etatrig",  "Lepton $\eta$", [0, 1.479, 2.4])),

########## Hists I added
        'ht_atlas'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("ht_atlas",      "H$_{T}$ ATLAS (GeV)", 24, 5, 365)),
        'mtw'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("mtw",      "m$_{T}^{W}$ (GeV)", 10, 0, 200)),
        'met_mtw'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("met_mtw",      "E$_{T}^{miss}+$m$_{T}^{W}$ (GeV)", 6, 70, 200)),
        'mlb'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("mlb",      "m$_{\ell b}$ (GeV)", 11, 5, 170)),
        'beta'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("beta",      "B-jet $\eta$ ", 8, -2.5, 2.50)),
        'deltaeta'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("deltaeta",      "Jets $\Delta \eta$ ", 8, 0, 2.50)),
        'topmass'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("topmass",      " $M_{top}$ ", 20, 100, 500)),        

        'mtwnocut'   : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("mtwnocut",     "m$_{T}^{W}$ (GeV)", 6, 0, 120)),
        'ht_atlas_nocut'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("ht_atlas_nocut",      "H$_{T}$ ATLAS (GeV)", 10, 0, 340)),


#New ones used for training (according to references)
        'u0mass'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("u0mass",      "$m_{u}$ (GeV)", 10, 0, 30)),
        'abseeta'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("abseeta",      "Electron $|\eta|$ ", 8, 0, 2.5)),
        'absmeta'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("absmeta",      "Muon $|\eta|$ ", 8, 0, 2.5)),
        'drub'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("drub",      " $\Delta R_{u,b}$ ", 14, 0, 6)),
#        'coste'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("coste",      " cos($\Theta^{*}$) ", 15, -1, 1)), 
#        'coste_top'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("coste_top",      " cos($\Theta^{*}$) ", 15, -1, 1)),  
#        'costm'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("costm",      " cos($\Theta^{*}$) ", 15, -1, 1)), 
#        'costm_top'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("costm_top",      " cos($\Theta^{*}$) ", 15, -1, 1)), 
        
#        'coste_1'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("coste_1",      " cos($\Theta^{*}$) ", 15, -1, 1)), 
        'coste_2'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("coste_2",      " cos($\Theta^{*}$) ", 15, -1, 1)), 
#        'coste_3'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("coste_3",      " cos($\Theta^{*}$) ", 15, -1, 1)),                         

#        'costm_1'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("costm_1",      " cos($\Theta^{*}$) ", 15, -1, 1)), 
        'costm_2'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("costm_2",      " cos($\Theta^{*}$) ", 15, -1, 1)), 
#        'costm_3'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("costm_3",      " cos($\Theta^{*}$) ", 15, -1, 1)),          
        
        'absweta'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("absweta",      " W boson $|\eta|$ ", 20, 0, 3.5)), 
        'topeta'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("topeta",      " top quark $\eta$ ", 8, -2.5, 2.5)), 
       
        'mub'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("mub",     "m(u,b) (GeV)", 10, 0, 200)),
        'DptWub'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("DptWub",     "|$\Delta_{p_{T}}$(W,ub)| (GeV)", 10, 0, 200)),
        'DphiWub'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("DphiWub",     "|$\Delta\phi$(W,ub)| ", 10, 0, 5)),
        'Detalu'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("Detalu",     "|$\Delta\eta(\ell,u)$| ", 10, 0, 5)),
#        'Detalb'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("Detalb",     "|$\Delta\eta(\ell,b)$| ", 10, 0, 5)),
#        'DRlu'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("DRlu",     "|$\Delta R(\ell,u)$| ", 10, 0, 5)),
########## hists i added end




        # processor.column_accumulator
        # Regions: 
        #  A: 2j1b, 3j1b, 3j2b,   
        #  B: 3j1b
        # Vars:
        # njets, nbtags, ht, dRlb, st, ptlb, sumallpt, ptuu, muu, medianDRuu, u0pt, u0eta, ptjj, mjj, medianDRjj, minDRjj, mlb, mt
        '2j1b_njets' : processor.column_accumulator(np.array([])),
        '2j1b_nbtags' : processor.column_accumulator(np.array([])),
        '2j1b_ht' : processor.column_accumulator(np.array([])),
        '2j1b_st' : processor.column_accumulator(np.array([])),
        '2j1b_sumAllPt' : processor.column_accumulator(np.array([])),
        '2j1b_leta' : processor.column_accumulator(np.array([])), #ojo
        '2j1b_lpt' : processor.column_accumulator(np.array([])), #ojo
        '2j1b_j0pt' : processor.column_accumulator(np.array([])),
        '2j1b_j0eta' : processor.column_accumulator(np.array([])),
        '2j1b_u0pt' : processor.column_accumulator(np.array([])),
        '2j1b_u0eta' : processor.column_accumulator(np.array([])),
        '2j1b_ptjj' : processor.column_accumulator(np.array([])),
        '2j1b_mjj' : processor.column_accumulator(np.array([])),
        '2j1b_medianDRjj' : processor.column_accumulator(np.array([])),
        '2j1b_minDRjj' : processor.column_accumulator(np.array([])),
        '2j1b_mlb' : processor.column_accumulator(np.array([])),
        '2j1b_mt' : processor.column_accumulator(np.array([])),
        '2j1b_ptsumveclb' : processor.column_accumulator(np.array([])),
        '2j1b_drlb' : processor.column_accumulator(np.array([])),
        '2j1b_druu' : processor.column_accumulator(np.array([])),
        '2j1b_druumedian' : processor.column_accumulator(np.array([])),
        '2j1b_muu' : processor.column_accumulator(np.array([])),
        '2j1b_ptuu' : processor.column_accumulator(np.array([])),
        
        ##Added for tchannel
        '2j1b_absu0eta' : processor.column_accumulator(np.array([])),
        '2j1b_ht_atlas' : processor.column_accumulator(np.array([])),
        '2j1b_mtw' : processor.column_accumulator(np.array([])),
        '2j1b_beta' : processor.column_accumulator(np.array([])),
        '2j1b_deltaeta' : processor.column_accumulator(np.array([])),
        '2j1b_topmass' : processor.column_accumulator(np.array([])),
        '2j1b_u0mass' : processor.column_accumulator(np.array([])),
        '2j1b_absleta' : processor.column_accumulator(np.array([])), #ojo con esta
        '2j1b_topeta' : processor.column_accumulator(np.array([])),
        '2j1b_absweta' : processor.column_accumulator(np.array([])),
        '2j1b_drub' : processor.column_accumulator(np.array([])),
        '2j1b_cost2' : processor.column_accumulator(np.array([])),   #ojo con esta tb
        
        '2j1b_mub' : processor.column_accumulator(np.array([])),
        '2j1b_DptWub' : processor.column_accumulator(np.array([])),
        '2j1b_DphiWub' : processor.column_accumulator(np.array([])),
        '2j1b_Detalu' : processor.column_accumulator(np.array([])),
#        '2j1b_Detalb' : processor.column_accumulator(np.array([])),
#        '2j1b_DRlu' : processor.column_accumulator(np.array([])),
        


        })

    @property
    def accumulator(self):
        return self._accumulator

    @property
    def columns(self):
        return self._columns

    # Main function: run on a given dataset
    def process(self, events):
        # Dataset parameters
        dataset = events.metadata["dataset"]
        histAxisName = self._samples[dataset]["histAxisName"]
        histAxisName_control = self._samples[dataset]["histAxisName"]
        year         = self._samples[dataset]["year"]
        xsec         = self._samples[dataset]["xsec"]
        sow          = self._samples[dataset]["nSumOfWeights"]
        isData       = self._samples[dataset]["isData"]
        isSystSample = ('mtop' in histAxisName) or ('hdamp' in histAxisName) or ('UE' in histAxisName)
        doPS         = (histAxisName in ['tchan', 'tchannel','tbarchannel']) and events.PSWeight is not None and len(events.PSWeight[0])>=4
        doPDFunc = "sumPDFWeights" in self._samples[dataset]

        # Get the lumi mask for 5 TeV data
        golden_json_path = cafea_path("data/goldenJsons/Cert_306546-306826_5TeV_EOY2017ReReco_Collisions17_JSON.txt")

        if doPDFunc:
          sowPDF       = self._samples[dataset]["sumPDFWeights"]
          sowScale     = self._samples[dataset]["sumScaleWeights"]
          PDFnorm = 1./np.array(sowPDF)
          Scalenorm = 1./np.array(sowScale)
          scaleweights      = events.LHEScaleWeight.to_numpy()
          scaleweights_bins = ak.local_index(events.LHEScaleWeight)
          pdfweights        = events.LHEPdfWeight.to_numpy()
          pdfweights_bins   = ak.local_index(events.LHEPdfWeight)
          scaleweights      = scaleweights * Scalenorm
          pdfweights        = pdfweights * PDFnorm



        # Initialize objects
        met  = events.MET
        e    = events.Electron
        mu   = events.Muon
        tau  = events.Tau
        jets = events.Jet
        
        
        

        # Pre-selection (must be updated with 5TeV definitions)
        #e["idEmu"] = ttH_idEmu_cuts_E3(e.hoe, e.eta, e.deltaEtaSC, e.eInvMinusPInv, e.sieie)
        #e["conept"] = coneptElec(e.pt, e.mvaTTH, e.jetRelIso)
        #mu["conept"] = coneptMuon(mu.pt, mu.mvaTTH, mu.jetRelIso, mu.mediumId)
        e["btagDeepB"] = ak.fill_none(e.matched_jet.btagDeepB, -99)
        mu["btagDeepB"] = ak.fill_none(mu.matched_jet.btagDeepB, -99)

        # Muon selection
        mu["isLoose"] = MuonLoose(mu.pt, mu.eta, mu.dxy, mu.dz, mu.sip3d, mu.mediumPromptId, mu.btagDeepB, ptCut=20, etaCut=2.4)
        mu["isMVA"]= MuonMVA(mu.miniPFRelIso_all, mu.mvaTTH)
        #mu["isIso"] = MuonTightIso(mu)
        #mu["isTight"] = isMuonPOGMnoIso(mu, ptCut=20, etaCut=2.4) # This is medium prompt, includes dxy, dz

        # Electron selection
        #e['pt'] = GetElecPt(e.pt, e.eta, isdata = isData)
        #e['pt'] = GetElecPtSmear(e.pt, e.eta, isdata = isData)
        GetElecScale5TeV(e, run=306936, isData=isData)
        e['isLoose'] = ElecLoose(e.pt, e.eta, e.lostHits, e.sip3d, e.dxy, e.dz, e.btagDeepB, e.convVeto, e.mvaFall17V2noIso_WPL, 20, 2.4)
        e['isMVA']   = ElecMVA(e.miniPFRelIso_all, e.mvaTTH)
        #e['isTight'] = (e.mvaFall17V2Iso_WP80)&(e.pt>20)&(np.abs(e.eta)<2.4)
        #e['isNoIso'] = (e.mvaFall17V2noIso_WP80)&(e.pt>20)&(np.abs(e.eta)<2.4)&(e.isTight == 0)
        

        # Build loose collections
        m_sel = mu[mu.isLoose & mu.isMVA]
        e_sel = e[e.isLoose & e.isMVA]
        m_fake = mu[mu.isLoose & (mu.isMVA == 0)]
        e_fake = e[e.isLoose & (e.isMVA == 0)]
        #m_sel = mu[mu.isTight & mu.isIso]
        #m_fake = mu[mu.isTight & (mu.isIso == 0)]
        #e_sel = e[e.isTight]
        #e_fake = e[e.isNoIso]
        e0 = e_sel[ak.argmax(e_sel.pt,axis=-1,keepdims=True)]
        m0 = m_sel[ak.argmax(m_sel.pt,axis=-1,keepdims=True)]
        
      
        
        e_fake0 = e_fake[ak.argmax(e_fake.pt,axis=-1,keepdims=True)]
        m_fake0 = m_fake[ak.argmax(m_fake.pt,axis=-1,keepdims=True)]
  
        if not isData:
          AttachElectronSF(e_sel,year='5TeV')
          AttachMuonSF(m_sel,year='5TeV')
          #AttachTrigSF(e0, m0, events)

        l_sel = ak.with_name(ak.concatenate([e_sel, m_sel], axis=1), 'PtEtaPhiMCandidate')
        leps  = ak.with_name(ak.concatenate([e_sel, m_sel], axis=1), 'PtEtaPhiMLorentzVector')
        fakes = ak.with_name(ak.concatenate([e_fake, m_fake], axis=1), 'PtEtaPhiMLorentzVector')
        l_sel_padded = ak.pad_none(l_sel, 1)
        lsel0 = l_sel_padded[:,0]
        

        l_fake = ak.with_name(ak.concatenate([e_fake, m_fake], axis=1), 'PtEtaPhiMCandidate')
        l_fake_padded = ak.pad_none(l_fake, 1)
        lfake0 = l_fake_padded[:,0]
        
        
        
        events.MET['pt_raw'] = events.RawMET.pt

        if not isData:
          e_sel_padded = ak.pad_none(e_sel, 1)
          m_sel_padded = ak.pad_none(m_sel, 1)
          events['sf_e'] = ak.fill_none(e_sel_padded[:,0].sf_nom, 1)
          events['sf_e_hi'] = ak.fill_none(e_sel_padded[:,0].sf_hi, 1)
          events['sf_e_lo'] = ak.fill_none(e_sel_padded[:,0].sf_lo, 1)
          events['sf_m'] = ak.fill_none(m_sel_padded[:,0].sf_nom, 1)
          events['sf_m_hi'] = ak.fill_none(m_sel_padded[:,0].sf_hi, 1)
          events['sf_m_lo'] = ak.fill_none(m_sel_padded[:,0].sf_lo, 1)
          AddSFs(events, l_sel) # for 2l

        events['isem'] = (ak.num(m_sel) == 1) & (ak.num(e_sel) == 1)
        events['ismm'] = (ak.num(m_sel) == 2) & (ak.num(e_sel) == 0)
        events['isee'] = (ak.num(m_sel) == 0) & (ak.num(e_sel) == 2)
        events['ise' ] = (ak.num(m_sel) == 0) & (ak.num(e_sel) == 1)
        events['ism' ] = (ak.num(m_sel) == 1) & (ak.num(e_sel) == 0)
        events['ise_fake' ] = (ak.num(m_sel) == 0) & ((ak.num(e_sel)) == 0) & (ak.num(e_fake) == 1)
        events['ism_fake' ] = (ak.num(m_sel) == 0) & ((ak.num(e_sel)) == 0) & (ak.num(m_fake) == 1)
        
        
        #events['isl_plus']=((ak.num(m_sel) == 0) & (ak.num(e_sel) == 1)&(e_sel.pdgId>0))|((ak.num(m_sel) == 1) & (ak.num(e_sel) == 0)&(m_sel.pdgId>0))
        
        
        
        if not isData:
          AttachTrigSF(e0, m0, events)
        # Jet cleaning, before any jet selection
        vetos_tocleanjets = ak.with_name( l_sel, "PtEtaPhiMCandidate")
        tmp = ak.cartesian([ak.local_index(jets.pt), vetos_tocleanjets.jetIdx], nested=True)
        cleanedJets = jets[~ak.any(tmp.slot0 == tmp.slot1, axis=-1)] # this line should go before *any selection*, otherwise lep.jetIdx is not aligned with the jet index

        if not isData:
          cleanedJets["pt_raw"] = (1 - cleanedJets.rawFactor)*cleanedJets.pt
          cleanedJets["mass_raw"] = (1 - cleanedJets.rawFactor)*cleanedJets.mass
          cleanedJets["pt_gen"] = ak.values_astype(ak.fill_none(cleanedJets.matched_gen.pt, 0), np.float32)
          cleanedJets["rho"] = ak.broadcast_arrays(events.fixedGridRhoFastjetAll, cleanedJets.pt)[0]
          events_cache = events.caches[0]
          corrected_jets = cleanedJets
          corrected_jets = jet_factory.build(cleanedJets, lazy_cache=events_cache)
          cleanedJets = corrected_jets
          met = met_factory.build(events.MET, corrected_jets, events.caches[0])

        else:
          cleanedJets["pt_raw"] = (1 - cleanedJets.rawFactor)*cleanedJets.pt
          cleanedJets["mass_raw"] = (1 - cleanedJets.rawFactor)*cleanedJets.mass
          cleanedJets["rho"] = ak.broadcast_arrays(events.fixedGridRhoFastjetAll, cleanedJets.pt)[0]
          corrected_jets = jet_factory_data.build(cleanedJets, lazy_cache=events.caches[0])
          cleanedJets = corrected_jets
          events_cache = events.caches[0]
          met = met_factory.build(events.MET, corrected_jets, events_cache)

        ################################ All this depends on jet pt
        jetptcut = 25
        jetetacut=2.4
        
        # Recommendations for 94X: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
        # Medium DeepCSV for 94X wp = 0.4941 , tight = 0.8001, loose = 0.1522
        # Medium DeepJet for 94X wp = 0.3033 
        # We're using DeepCSV...
        wp = 0.4941#0.8001 #0.4941 # medium
        goodJets, metpt, btagSF = GetJESMETvar(corrected_jets, met, jetptcut, btagwp=wp, isData=isData,var='')
        goodJets_norm = goodJets
        metpt_norm = metpt
        
      
        ######### SFs, weights, systematics ##########
        # Btag SF following 1a) in https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods 
        if not isData:
          btagSF, btagSFUp, btagSFDo = GetBtagSF5TeV(goodJets.pt, goodJets.eta, goodJets.hadronFlavour, goodJets.isBtag, True) #


        ### Trigger
        
        trigem = (events.HLT.HIMu17) | (events.HLT.HIEle15_WPLoose_Gsf)
        trigee = (events.HLT.HIEle15_WPLoose_Gsf) | (events.HLT.HIEle17_WPLoose_Gsf)
        trigmm = (events.HLT.HIMu17)
        trige  = (events.HLT.HIEle20_WPLoose_Gsf)
        trigm  = events.HLT.HIMu17#HIL3Mu20
        passtrige = (events.HLT.HIEle20_WPLoose_Gsf)
        passtrigm = events.HLT.HIMu17#HIL3Mu20
        '''
        #tchannel
        trigem = np.ones_like(events['event'], dtype=bool)#(events.HLT.HIMu17) | (events.HLT.HIEle15_WPLoose_Gsf)
        trigee = np.ones_like(events['event'], dtype=bool)#(events.HLT.HIEle15_WPLoose_Gsf) | (events.HLT.HIEle17_WPLoose_Gsf)
        trigmm = np.ones_like(events['event'], dtype=bool)#(events.HLT.HIMu17)
        trige  = np.ones_like(events['event'], dtype=bool)#(events.HLT.HIEle20_WPLoose_Gsf)
        trigm  = np.ones_like(events['event'], dtype=bool)#events.HLT.HIMu17#HIL3Mu20
        passtrige = np.ones_like(events['event'], dtype=bool)#(events.HLT.HIEle20_WPLoose_Gsf)
        passtrigm = np.ones_like(events['event'], dtype=bool)#events.HLT.HIMu17#HIL3Mu20
        '''

        # Single electron events: trige, only from HighEGJet
        # Single muon events: trigm, only from SingleMuon
        # ee events: trigee, only from HighEGJet
        # mm events: trigmm, only from SingleMuon
        # em events: in SingleMuon: pass trigmm, in HighEGJet: pass trigee and not trigmm

        if isData:
          if   histAxisName=='HighEGJet': 
            trigem = ((events.HLT.HIEle15_WPLoose_Gsf) | (events.HLT.HIEle17_WPLoose_Gsf)) & ((events.HLT.HIMu17)==0)
            trigee = (events.HLT.HIEle15_WPLoose_Gsf) | (events.HLT.HIEle17_WPLoose_Gsf)
            trigmm = np.zeros_like(events['event'], dtype=bool)
            trige  = events.HLT.HIEle20_WPLoose_Gsf
            trigm  = np.zeros_like(events['event'], dtype=bool)
          elif histAxisName=='SingleMuon': 
            trigem = events.HLT.HIMu17
            trigee = np.zeros_like(events['event'], dtype=bool)
            trigmm = (events.HLT.HIMu17)
            trige  = np.zeros_like(events['event'], dtype=bool)
            trigm  = events.HLT.HIMu17#HIL3Mu20

        # We need weights for: normalization, lepSF, triggerSF, pileup, btagSF...
        # prefire
        if not isData:
          from cafea.modules.prefire import GetPrefireWeights
          # Must be done before any selection!
          prefweight     = GetPrefireWeights(events.Jet, events.Photon, events.Electron, var= 0)
          prefweightUp   = GetPrefireWeights(events.Jet, events.Photon, events.Electron, var= 1)
          prefweightDown = GetPrefireWeights(events.Jet, events.Photon, events.Electron, var=-1)

        weights_dict = {}
        if (isData): genw = np.ones_like(events["event"])
        else:        genw = events["genWeight"]
        for ch_name in ["em", "e", "m", 'ee', 'mm']:
          weights_dict[ch_name] = coffea.analysis_tools.Weights(len(events),storeIndividual=True)
          weights_dict[ch_name].add("norm",genw if isData else (xsec/sow)*genw)
          
          if not isData: # Apply SFs
            if ch_name in ["em", "ee", "mm"]:
              weights_dict[ch_name].add("lepSF", events.sf_2l, events.sf_2l_hi, events.sf_2l_lo)
            else:
              #weights_dict[ch_name].add("lepSF", ak.copy(events.sf_1l), ak.copy(events.sf_1l_hi), ak.copy(events.sf_1l_lo))
              weights_dict[ch_name].add("elecSF", ak.copy(events.sf_e), ak.copy(events.sf_e_hi), ak.copy(events.sf_e_lo))
              weights_dict[ch_name].add("muonSF", ak.copy(events.sf_m), ak.copy(events.sf_m_hi), ak.copy(events.sf_m_lo))
            weights_dict[ch_name].add("trigSF", ak.copy(events.sf_trig), ak.copy(events.sf_trig_hi), ak.copy(events.sf_trig_lo))
            weights_dict[ch_name].add("btagSF", ak.copy(btagSF), ak.copy(btagSFUp), ak.copy(btagSFDo))
            weights_dict[ch_name].add("prefire", ak.copy(prefweight), ak.copy(prefweightUp), ak.copy(prefweightDown))
        
     
          # PS = ISR, FSR (on tchan only, see line 410)
          if doPS: 
            i_ISRdown = 0; i_FSRdown = 1; i_ISRup = 2; i_FSRup = 3
            ISRUp = ak.flatten(events.PSWeight[ak.Array(ak.local_index(events.PSWeight)==i_ISRup)])
            ISRDo = ak.flatten(events.PSWeight[ak.Array(ak.local_index(events.PSWeight)==i_ISRdown)])
            FSRUp = ak.flatten(events.PSWeight[ak.Array(ak.local_index(events.PSWeight)==i_FSRup)])
            FSRDo = ak.flatten(events.PSWeight[ak.Array(ak.local_index(events.PSWeight)==i_FSRdown)])
            weights_dict[ch_name].add('ISR', np.ones_like(events["event"]), ISRUp, ISRDo)
            weights_dict[ch_name].add('FSR', np.ones_like(events["event"]), FSRUp, FSRDo)

        # Add systematics
        systList = ["norm"]
        systJES = ['JER','MC', 'AbsStat', 'AbsScale', 'AbsMPF', 'Frag', 'ECAL', 'HCAL', 'Flavor', 'RelStat', 'RelPt', 'RelBal', 'RelJER', 'L3Res','MET_UnclusteredEnergy','Total'] if splitJES else ['Total','JER','MET_UnclusteredEnergy']
        systJets = [x + 'Up' for x in systJES] + [x + 'Down' for x in systJES]
        if not isData and not isSystSample: systList = systList + ["elecSFUp","elecSFDown", "muonSFUp", "muonSFDown", "btagSFUp", "btagSFDown", "prefireUp", "prefireDown", "trigSFUp", "trigSFDown"]+systJets#, "trigSFUp", "trigSFDown"] + systJets
        if doPS: systList += ['ISRUp', 'ISRDown', 'FSRUp', 'FSRDown']
        if not doSyst or isData or isSystSample: systList = ["norm"]

        # Add selections...
        selections = PackedSelection(dtype='uint64')
        selections.add("em", ( (events.isem)&(trigem)))
        selections.add("ee", ( (events.isee)&(trigee)))
        selections.add("mm", ( (events.ismm)&(trigmm)))
        selections.add("e", ( (events.ise)&(trige)))
        selections.add("m", ( (events.ism)&(trigm)))
        selections.add("e_fake", ( (events.ise_fake)&(trige)))
        selections.add("m_fake", ( (events.ism_fake)&(trigm)))
        
       # selections.add('lep_plus',((events.isl_plus)&((trige)|(trigm))))

        selections.add("incl", ak.ones_like(met.pt, dtype=bool))
        metfilters = PassMETfilters(events,isData)
        selections.add("METfilters", metfilters)
        # Counts
        counts = np.ones_like(events['event'], dtype=float)
 
        # Initialize the out object
        hout = self.accumulator.identity()
        #channels = ['e', 'm', 'e_fake', 'm_fake'] #['em', 'e', 'm', 'ee', 'mm', 'e_fake', 'm_fake'] 
        channels=['e','m','e_fake','m_fake']#,'mm']
        
        #levels = ['incl', 'g2jets', 'g3jets', 'g4jets', '0b','3j0b','4j0b', 'g5j0b', '2j1b', '3j1b', '3j2b','4j1b', '4j2b', 'g5j1b', 'g5j2b']
        levels = ['2j1b','3j1b','3j2b']

        # For trigger studies
        e0pt      = ak.flatten(e0[ak.num(e_sel)>=1].pt)
        e0pt_pass = ak.flatten(e0[(ak.num(e_sel)>=1)&(trige)].pt)
        m0pt     = ak.flatten(m0[ak.num(m_sel)>=1].pt)
        m0pt_pass = ak.flatten(m0[(ak.num(m_sel)>=1)&(trigm)].pt)
        hout['pttrig'].fill(sample=histAxisName, channel='e', level='allden', pttrig=e0pt, weight=np.ones_like(e0pt))
        hout['pttrig'].fill(sample=histAxisName, channel='e', level='allnum', pttrig=e0pt_pass, weight=np.ones_like(e0pt_pass))
        hout['pttrig'].fill(sample=histAxisName, channel='m', level='allden', pttrig=m0pt, weight=np.ones_like(m0pt))
        hout['pttrig'].fill(sample=histAxisName, channel='m', level='allnum', pttrig=m0pt_pass, weight=np.ones_like(m0pt_pass))
        # eµ cross-trigger
        trem_e0pt       = ak.flatten(e0[(events.isem)&(passtrigm)].pt)        
        trem_e0pt_pass  = ak.flatten(e0[(events.isem)&(passtrigm)&(passtrige)].pt)
        trem_m0pt       = ak.flatten(m0[(events.isem)&(passtrige)].pt)
        trem_m0pt_pass  = ak.flatten(m0[(events.isem)&(passtrigm)&(passtrige)].pt)
        trem_e0eta      = np.abs(ak.flatten(e0[(events.isem)&(passtrigm)].eta))
        trem_e0eta_pass = np.abs(ak.flatten(e0[(events.isem)&(passtrigm)&(passtrige)].eta))
        trem_m0eta      = np.abs(ak.flatten(m0[(events.isem)&(passtrige)].eta))
        trem_m0eta_pass = np.abs(ak.flatten(m0[(events.isem)&(passtrigm)&(passtrige)].eta))
        hout['pttrig'].fill(sample=histAxisName, channel='e', level='den', pttrig=trem_e0pt, weight=np.ones_like(trem_e0pt))
        hout['pttrig'].fill(sample=histAxisName, channel='e', level='num', pttrig=trem_e0pt_pass, weight=np.ones_like(trem_e0pt_pass))
        hout['pttrig'].fill(sample=histAxisName, channel='m', level='den', pttrig=trem_m0pt, weight=np.ones_like(trem_m0pt))
        hout['pttrig'].fill(sample=histAxisName, channel='m', level='num', pttrig=trem_m0pt_pass, weight=np.ones_like(trem_m0pt_pass))
        hout['etatrig'].fill(sample=histAxisName, channel='e', level='den', etatrig=trem_e0eta, weight=np.ones_like(trem_e0eta))
        hout['etatrig'].fill(sample=histAxisName, channel='e', level='num', etatrig=trem_e0eta_pass, weight=np.ones_like(trem_e0eta_pass))
        hout['etatrig'].fill(sample=histAxisName, channel='m', level='den', etatrig=trem_m0eta, weight=np.ones_like(trem_m0eta))
        hout['etatrig'].fill(sample=histAxisName, channel='m', level='num', etatrig=trem_m0eta_pass, weight=np.ones_like(trem_m0eta_pass))
        

        # These are computed once for all the non-JES systematics
        njets_nom, nbtags_nom, nujets_nom = GetNjetNbtagsNujets(goodJets)
        nbtagnjets = GetNBtagNJets(njets_nom, nbtags_nom)
        ht_nom = ak.sum(goodJets.pt,axis=-1)
        j0_nom, drjj_nom, drjjmedian_nom, mjj_nom, ptjj_nom = GetJetVariables(goodJets)
        b0_nom, drbb_nom, drbbmedian_nom, mbb_nom, ptbb_nom = GetJetVariables(goodJets[(goodJets.isBtag)])
        bhard_nom=GetMaxBJet(goodJets)
        u0_nom, druu_nom, druumedian_nom, muu_nom, ptuu_nom = GetJetVariables(goodJets[(goodJets.isBtag==0)])
        ptSumVecAll_nom, ptSumVeclb_nom, dRlb_nom, st_nom   = GetJetLepVar(goodJets, leps)
        ptSumVecAll_fak, ptSumVeclb_fak, dRlb_fak, st_fak   = GetJetLepVar(goodJets, fakes)

        


        # Loop over the hists we want to fill
        for syst in systList:
          j0, drjj, drjjmedian, mjj, ptjj   = (j0_nom, drjj_nom, drjjmedian_nom, mjj_nom, ptjj_nom)
          b0, drbb, drbbmedian, mbb, ptbb   = (b0_nom, drbb_nom, drbbmedian_nom, mbb_nom, ptbb_nom)
          bhard=bhard_nom
          u0, druu, druumedian, muu, ptuu   = (u0_nom, druu_nom, druumedian_nom, muu_nom, ptuu_nom)
          ptSumVecAll, ptSumVeclb, dRlb, st = (ptSumVecAll_nom, ptSumVeclb_nom, dRlb_nom, st_nom)
          njets, nbtags, nujets = (njets_nom, nbtags_nom, nujets_nom)
          ht_var = ht_nom
          met_pt = metpt
          njets_var = njets
          nbtags_var = nbtags
          nujets_var = nujets
          goodJets = goodJets_norm
          metpt = metpt_norm

          if syst in systJets:
            goodJets, metpt, btagSF = GetJESMETvar(corrected_jets, met, jetptcut, btagwp=wp, isData=isData, var=syst)
            njets, nbtags, nujets = GetNjetNbtagsNujets(goodJets)
            njets_var = njets
            nbtags_var = nbtags
            nujets_var = nujets
            ht_var = ak.sum(goodJets.pt,axis=-1)
            met_pt = metpt
            j0, drjj, drjjmedian, mjj, ptjj   = GetJetVariables(goodJets)
            b0, drbb, drbbmedian, mbb, ptbb   = GetJetVariables(goodJets[(goodJets.isBtag)])
            bhard=GetMaxBJet(goodJets)
            u0, druu, druumedian, muu, ptuu   = GetJetVariables(goodJets[(goodJets.isBtag==0)])
            ptSumVecAll, ptSumVeclb, dRlb, st = GetJetLepVar(goodJets, leps)

          jet0pt  = ak.flatten(j0.pt)
          jet0eta = ak.flatten(j0.eta)
         
          
          
          b0pt  = ak.flatten(b0.pt)
          b0eta = ak.flatten(b0.eta)
          b0phi=ak.flatten(b0.phi)
          
          
          u0pt  = ak.flatten(u0.pt)
          u0eta = ak.flatten(u0.eta)
          u0phi=ak.flatten(u0.phi)
          u0mass=ak.flatten(u0.mass)
          
          
          dRub=np.sqrt((u0eta-b0eta)**2+(u0phi-b0phi)**2)

          for ch in channels:
            if syst in ['elecSFUp', 'elecSFDown', 'muonSFUp', 'muonSFDown'] and ch in ['ee', 'mm', 'em']: continue
            if ch in ['e_fake', 'm_fake']:
              ptSumVecAll, ptSumVeclb, dRlb, st = (ptSumVecAll_fak, ptSumVeclb_fak, dRlb_fak, st_fak)
            for lev in levels:
              cutschan  = [ch] + ['incl']+ ['METfilters']
              cutsel = selections.all(*cutschan)
              
              var=ch+'0'
              
              
				  
                
#######################################################################################################    Adding by myself extra cuts to replicate ATLAS selection
#              cutbetal=abs(b0eta)<2.5
 #             cutbeta=ak.Array(cutbetal)
  #            cutbeta=ak.fill_none(cutbeta,False)    #For adding the eta cut in b jet (me)
   #           
    #          maskcutuetal = (1.5 < np.abs(u0eta)) & (np.abs(u0eta) < 4.0)
     #         cutueta=ak.fill_none(maskcutuetal,False)    #For adding the eta cut in non-b jet (me)
              
              
#              cutdeltaetal=abs(b0eta-u0eta)>1.5
 #             cutdeltaeta=ak.Array(cutdeltaetal)
  #            cutdeltaeta=ak.fill_none(cutdeltaeta,False)  #For adding the difference in eta bt b and u jets cut (me)
            
              
             
              lep_pt=getattr(eval(var),'pt')
              #lep_px=getattr(eval(var),'px')
              #lep_py=getattr(eval(var),'py')
              #lep_pz=getattr(eval(var),'pz')
              lep_energy=getattr(eval(var),'energy')
              lep_phi=getattr(eval(var),'phi')
              lep_eta=getattr(eval(var),'eta')
              lep_mass=getattr(eval(var),'mass')
              
          
              
              ht_atlas=ht_var+lep_pt+met_pt
              

              
              mtw=np.sqrt(2*lep_pt*met_pt*(1-np.cos(abs(lep_phi-met.phi))))   #variable atlas defines in their paper
              
              mtw=ak.fill_none(mtw,False)
              mtw = ak.to_numpy(ak.flatten(mtw))
             
              cost=np.cos((2*np.arctan(np.exp(-u0eta)))-(2*np.arctan(np.exp(-ak.flatten(lep_eta)))))

              

            
              
              cutmtwl=mtw>50
              cutmtw=ak.Array(cutmtwl)   #For adding the cut in Mw (me)
              cutmtw=ak.fill_none(cutmtw,False)
              
              
#              cutmet_mtw=met_pt+mtw>70  #For adding cut in sum of MET and Mw (me)
 #             cutmet_mtw=ak.fill_none(cutmet_mtw,False)
              
              cutht=ht_atlas>170  #For adding cut in Ht (me) Beware this is a different definition from our usual Ht
              cutht=ak.fill_none(cutht,False)
              
              
              if ch == 'e': mlb=GetMlb1b(e0,b0)
              if ch== 'm': mlb=GetMlb1b(m0,b0)
              if ch == 'e_fake': mlb=GetMlb1b(e_fake0,b0)
              if ch== 'm_fake': mlb=GetMlb1b(m_fake0,b0)
              
              cutmlb=mlb<125
              cutmlb=ak.fill_none(cutmlb,False)
              
              
              
              
     #### NOTE: in this selection I'm putting together CMS criteria for good objects (leptons and jets basically, and at least in tt) with the selection described by ATLAS for tchannel.
     ### the good-object criteria of ATLAS is not replicated but seems similar to CMS one. If this was not the case, and we wanted to 
     ### to replicate the full ATLAS, we should revisit the good-objects that we have chosen to be so       
              
                            
########################################################################################################                   Up to here the cut thing
                                         
              cuts      = [ch] + [lev] + ['metg30'] 
              cutjets = GetCutJets(cuts, syst, met_pt, njets_var, nbtags_var)  #met_pt en vez de mtw si queremos de verdad estimar qcd con met y no con mtw
              
              cut = (cutsel) & (cutjets) & (cutmtw) & (cutht) & (cutmlb) 

              cut=ak.flatten(cut)
              cut = np.array(cut, dtype=bool)												#main cut
              
              
              
              cutsnoMET = [ch] + [lev] 
              cutjetsnomet = GetCutJets(cutsnoMET, syst, met_pt, njets_var, nbtags_var, nujets_var)
              cutnomet = (cutsel) & (cutjetsnomet)  & (cutmtw) & (cutht) & (cutmlb) 
              cutnomet=ak.flatten(cutnomet)
              cutnomet = np.array(cutnomet, dtype=bool)										#special cuts for met validations   OJO QUE SE LLAMA MET por consistencia con el resto del codigo pero es mtw
              
              
              

                       
 
               
#####################################       Time to do the W-boson recosntruction (not to raise problems in the minimization first we have to apply the cut) 		###################################
              w=[None]*len(cut)
              
              
              if ch in ['e','m']:
                for i in range(len(cut)):				  
                  if cutnomet[i]: w[i]=reco(lsel0[i],met[i],plot=False,plotpath='here/',n=0) 
              else:      				  
                for i in range(len(cut)):				  
                  if cutnomet[i]: w[i]=reco(lfake0[i],met[i],plot=False,plotpath='here/',n=0) 

             
             #### Claculation of variables of top quark and W boson
             
              
              w_candidates= EventCollection(w)
           
              # Calculate the masses for all events
              w_masses = w_candidates.M()

              invariant_masses = np.sqrt((w_candidates.E() + b0.energy)**2 - (w_candidates.Px() + b0.px)**2 - (w_candidates.Py() + b0.py)**2 - (w_candidates.Pz() + b0.pz)**2)

              if lev in ['3j2b']: invariant_masses = np.sqrt((w_candidates.E() + bhard.energy)**2 - (w_candidates.Px() + bhard.px)**2 - (w_candidates.Py() + bhard.py)**2 - (w_candidates.Pz() + bhard.pz)**2)

              invariant_masses=ak.flatten(invariant_masses)
              
              w_etas=w_candidates.eta()
              
              top_eta=(ak.flatten(lep_eta)+b0eta)*0.5
              
              
#              cost_top=np.cos((2*np.arctan(np.exp(-abs(u0eta-top_eta))))-(2*np.arctan(np.exp(-abs(ak.flatten(lep_eta)-top_eta)))))
              
#              cos1=np.cos(2*np.arctan(np.exp(-ak.flatten(lep_eta)))-2*np.arctan(np.exp(-w_etas)))
              cos2=np.cos(2*np.arctan(np.exp(-(ak.flatten(lep_eta)-w_etas)))-2*np.arctan(np.exp(-(u0eta-w_etas))))
#              cos2_abs=np.cos(2*np.arctan(np.exp(-abs(ak.flatten(lep_eta)-w_etas)))-2*np.arctan(np.exp(-abs(u0eta-w_etas))))
              
              
###Vars ATLAS 13 TeV:
              Wpt=np.sqrt(w_candidates.Px()*w_candidates.Px()+w_candidates.Py()*w_candidates.Py())
              ubpt=np.sqrt((b0.px+u0.px)**2+(b0.py+u0.py)**2)
              DptWub=ak.flatten(np.abs(Wpt-ubpt))
    
              mub=ak.flatten((u0+b0).mass)
              
              ubphi=ak.flatten(np.arctan2(u0.py+b0.py,u0.px+b0.px))
              Wphi=w_candidates.phi()
              DphiWub=np.abs(Wphi-ubphi)
              
              lep_eta_flat=ak.flatten(lep_eta)
              Detalu=np.abs(u0eta-lep_eta_flat)
              
#              Detalb=np.abs(b0eta-lep_eta_flat)
              
#              lep_phi_flat=ak.flatten(lep_phi)
#              DRlu=np.sqrt((u0eta-lep_eta_flat)**2+(u0phi-lep_phi_flat)**2)

######################################  Redefine cuts from here if wanting to add something on the w-top system ##############################################################
              
              cut0bs = [a or b for a, b in zip((np.abs(u0eta) < 2.0), (nbtags_var!=0))]

              cuttop=(invariant_masses>130)&(invariant_masses<190) #230
              cuttop=ak.fill_none(cuttop,False)
              cut=cut&cuttop    #& cut0bs
              cutnomet=cutnomet&cuttop  #& cut0bs

             

              
              if lev in ['2j0b']:												
                cutueta = np.abs(u0eta) < 2.0
                cutueta=ak.fill_none(cutueta,False)
                cut = (cutsel) & (cutjets) & (cutmtw)  & (cutht) & (cutueta)
                cut=ak.flatten(cut)
                cut = np.array(cut, dtype=bool)
                cutnomet = (cutsel) & (cutjetsnomet)  & (cutmtw)  & (cutht) & (cutueta)
                cutnomet=ak.flatten(cutnomet)
                cutnomet = np.array(cutnomet, dtype=bool)

            


             
              weights = weights_dict[ch if not 'fake' in ch else ch[0]].weight(syst if not syst in (['norm']+systJets) else None)
              
              ht_atlas=ak.fill_none(ht_atlas,False)
              ht_atlas=ak.flatten(ht_atlas)

              if histAxisName_control in ['W0JetsToLNu','W1JetsToLNu','W2JetsToLNu','W3JetsToLNu']:
              
               heavy_cond=goodJets.hadronFlavour>=4
               light_cond=(goodJets.hadronFlavour!=4)&(goodJets.hadronFlavour!=5)
               
               heavy_cut = np.array([any(subarray) for subarray in heavy_cond])
               light_cut = np.array([all(subarray) for subarray in light_cond])
               
               cut_heavyjets= cut & heavy_cut
               cut_lightjets= cut & light_cut
               cutnomet_heavyjets= cutnomet & heavy_cut
               cutnomet_lightjets= cutnomet & light_cut
               
               wjets_flavour_dict={'cut_heavyjets':cut_heavyjets,'cut_lightjets':cut_lightjets,
                'cutnomet_heavyjets':cutnomet_heavyjets,'cutnomet_lightjets':cutnomet_lightjets,
                'weights_control':weights,
                'heavyjets':histAxisName[:2]+'_heavyjets','lightjets':histAxisName[:2]+'_lightjets'}

              jetcategories=['dummy']
              if histAxisName_control in ['W0JetsToLNu','W1JetsToLNu','W2JetsToLNu','W3JetsToLNu']:jetcategories=['heavyjets','lightjets']
              for cat in jetcategories:
               if cat!='dummy': cutname='cut_'+cat; cut=wjets_flavour_dict[cutname];cutnometname='cutnomet_'+cat; cutnomet=wjets_flavour_dict[cutnometname]; weights=wjets_flavour_dict['weights_control']; histAxisName=wjets_flavour_dict[cat];
               if fillAcc and not isData and syst=='norm' and ch in ['e', 'm'] and lev in ['2j1b']:
       
                hout[f'{lev}_njets'] = processor.column_accumulator(njets[cut].to_numpy())
                hout[f'{lev}_nbtags'] = processor.column_accumulator(nbtags[cut].to_numpy())
                hout[f'{lev}_ht'] = processor.column_accumulator(ht_var[cut].to_numpy())
                hout[f'{lev}_st'] = processor.column_accumulator(st[cut].to_numpy())
                hout[f'{lev}_sumAllPt'] = processor.column_accumulator(ptSumVecAll[cut].to_numpy())
                hout[f'{lev}_leta'] = processor.column_accumulator(ak.flatten(l_sel.eta[cut]).to_numpy())
                hout[f'{lev}_j0pt'] = processor.column_accumulator(jet0pt[cut].to_numpy())
                hout[f'{lev}_j0eta'] = processor.column_accumulator(jet0eta[cut].to_numpy())
                hout[f'{lev}_u0pt'] = processor.column_accumulator(u0pt[cut].to_numpy())
                hout[f'{lev}_u0eta'] = processor.column_accumulator(u0eta[cut].to_numpy())

                hout[f'{lev}_ptjj'] = processor.column_accumulator(ak.flatten(ptjj[cut]).to_numpy())
                hout[f'{lev}_mjj'] = processor.column_accumulator(ak.flatten(mjj[cut]).to_numpy())
                hout[f'{lev}_medianDRjj'] = processor.column_accumulator(ak.flatten(drjjmedian[cut]).to_numpy())
                hout[f'{lev}_minDRjj'] = processor.column_accumulator(ak.flatten(drjj[cut]).to_numpy())
                hout[f'{lev}_mlb'] = processor.column_accumulator(mlb[cut].to_numpy())
                hout[f'{lev}_mt'] = processor.column_accumulator( ak.flatten(GetMT(l_sel, met)[cut]) .to_numpy())
                hout[f'{lev}_ptsumveclb'] = processor.column_accumulator( ptSumVeclb[cut].to_numpy())
                hout[f'{lev}_drlb'] = processor.column_accumulator( ak.flatten(dRlb[cut]).to_numpy())
                
                hout[f'{lev}_absu0eta'] = processor.column_accumulator( abs(u0eta[cut]).to_numpy())
                hout[f'{lev}_ht_atlas'] = processor.column_accumulator( ht_atlas[cut].to_numpy())
                hout[f'{lev}_mtw'] = processor.column_accumulator( mtw[cut])
                hout[f'{lev}_beta'] = processor.column_accumulator( b0eta[cut].to_numpy())
                hout[f'{lev}_deltaeta'] = processor.column_accumulator( abs(b0eta-u0eta)[cut].to_numpy())
                hout[f'{lev}_topmass'] = processor.column_accumulator( invariant_masses[cut].to_numpy())
                hout[f'{lev}_u0mass'] = processor.column_accumulator( u0mass[cut].to_numpy())
                hout[f'{lev}_absleta'] = processor.column_accumulator( abs(ak.flatten(lep_eta))[cut].to_numpy())
                hout[f'{lev}_topeta'] = processor.column_accumulator( top_eta[cut].to_numpy())
                hout[f'{lev}_absweta'] = processor.column_accumulator( abs(w_etas)[cut])
                hout[f'{lev}_drub'] = processor.column_accumulator( dRub[cut].to_numpy())
                hout[f'{lev}_cost2'] = processor.column_accumulator( cos2[cut].to_numpy())
                
                hout[f'{lev}_DptWub'] = processor.column_accumulator( DptWub[cut].to_numpy())
                hout[f'{lev}_mub'] = processor.column_accumulator( mub[cut].to_numpy())
                hout[f'{lev}_DphiWub'] = processor.column_accumulator( DphiWub[cut].to_numpy())
                hout[f'{lev}_Detalu'] = processor.column_accumulator( Detalu[cut].to_numpy())
                #hout[f'{lev}_Detalb'] = processor.column_accumulator( Detalb[cut].to_numpy())
                #hout[f'{lev}_DRlu'] = processor.column_accumulator( DRlu[cut].to_numpy())
                
                
                if lev == '3j1b':
                  hout[f'{lev}_druu'] = processor.column_accumulator( ak.flatten(druu[cut]).to_numpy())
                  hout[f'{lev}_druumedian'] = processor.column_accumulator( ak.flatten(druumedian[cut]).to_numpy())
                  hout[f'{lev}_muu'] = processor.column_accumulator( ak.flatten(muu[cut]).to_numpy())
                  hout[f'{lev}_ptuu'] = processor.column_accumulator( ak.flatten(ptuu[cut]).to_numpy())

               # Fill met norm histos
            
               cut_metg20 = (cutnomet) & GetCutJets(['metg20'], syst, met_pt, njets_var, nbtags_var)   # de aqui   para looser 2 tengo 25 20 15
               cut_metl20 = (cutnomet) & GetCutJets(['metl20'], syst, met_pt, njets_var, nbtags_var)
               cut_metg25 = (cutnomet) & GetCutJets(['metg25'], syst, met_pt, njets_var, nbtags_var)
               cut_metl25 = (cutnomet) & GetCutJets(['metl25'], syst, met_pt, njets_var, nbtags_var)
               cut_metl15 = (cutnomet) & GetCutJets(['metl15'], syst, met_pt, njets_var, nbtags_var)
               cut_metg30 = (cutnomet) & GetCutJets(['metg30'], syst, met_pt, njets_var, nbtags_var)
               cut_metg35 = (cutnomet) & GetCutJets(['metg35'], syst, met_pt, njets_var, nbtags_var)
               cut_metg15 = (cutnomet) & GetCutJets(['metg15'], syst, met_pt, njets_var, nbtags_var)   # a aqui  mtw-> met_pt si queremos hacer el calculo de qcd con met y no mtw
               cut_metl30 = (cutnomet) & GetCutJets(['metl30'], syst, met_pt, njets_var, nbtags_var)
             
               
               weights_metg20 = weights[cut_metg20]; weights_metl20 = weights[cut_metl20]
               weights_metg25 = weights[cut_metg25]; weights_metl25 = weights[cut_metl25]
               weights_metg15 = weights[cut_metg15]; weights_metl15 = weights[cut_metl15]
               weights_metg30 = weights[cut_metg30]; weights_metl30 = weights[cut_metl30]
               hout['counts_metg20'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut_metg20], syst=syst, weight=weights_metg20)
               hout['counts_metl20'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut_metl20], syst=syst, weight=weights_metl20)
               hout['counts_metg25'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut_metg25], syst=syst, weight=weights_metg25)
               hout['counts_metl25'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut_metl25], syst=syst, weight=weights_metl25)
               hout['counts_metg15'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut_metg15], syst=syst, weight=weights_metg15)
               hout['counts_metl15'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut_metl15], syst=syst, weight=weights_metl15)
               hout['counts_metg30'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut_metg30], syst=syst, weight=weights_metg30)
               hout['counts_metl30'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut_metl30], syst=syst, weight=weights_metl30)

              
               

               ### met dist with no cut
               hout['metnocut'].fill(sample=histAxisName, channel=ch, level=lev, metnocut=met.pt[cutnomet], syst=syst, weight=weights[cutnomet])
               
             #  hout['mtwnocut'].fill(sample=histAxisName, channel=ch, level=lev, mtwnocut=mtw[cutnomet], syst=syst, weight=weights[cutnomet])
               #hout['ht_atlas_nocut'].fill(sample=histAxisName, channel=ch, level=lev, ht_atlas_nocut=flat_ht_atlas_nocut, syst=syst, weight=weights[cut_noht])

   
               ### We need to have 2 jets in order to calculate dijet observables
               dijet_cut = (cut) & (GetCutJets(['g2jets'], syst, met_pt, njets_var, nbtags_var, nujets_var))
               dijuu_cut = (cut) & (GetCutJets(['g2ujets'], syst, met_pt, njets_var, nbtags_var, nujets_var))
               dijet_cut = np.array(dijet_cut, dtype=bool)
               dijuu_cut = np.array(dijuu_cut, dtype=bool)
               weights_dijet = weights[dijet_cut]
               weights_dijuu = weights[dijuu_cut]
               fdrjjmed = ak.flatten(drjjmedian[dijet_cut])
               fdrjjmin = ak.flatten(drjj[dijet_cut])
               hout['medianDRjj'].fill(sample=histAxisName, channel=ch, level=lev, medianDRjj=fdrjjmed, syst=syst, weight=weights_dijet)
               hout['minDRjj'].fill(sample=histAxisName, channel=ch, level=lev, minDRjj=fdrjjmin, syst=syst, weight=weights_dijet)
               if fillAll:
                 fmjj = ak.flatten(mjj[dijet_cut])
                 fptjj = ak.flatten(ptjj[dijet_cut])
                 fmedDRuu = ak.flatten(druumedian[dijuu_cut])
                 fminDRuu = ak.flatten(druu[dijuu_cut])
                 fmuu = ak.flatten(muu[dijuu_cut])
                 fptuu = ak.flatten(ptuu[dijuu_cut])
                 fu0pt = (u0pt[dijuu_cut])
                 fu0eta = (u0eta[dijuu_cut])
                 hout['mjj'].fill(sample=histAxisName, channel=ch, level=lev, mjj=fmjj, syst=syst, weight=weights_dijet)
                 hout['ptjj'].fill(sample=histAxisName, channel=ch, level=lev, ptjj=fptjj, syst=syst, weight=weights_dijet)
                 hout['medianDRuu'].fill(sample=histAxisName, channel=ch, level=lev, medianDRuu=fmedDRuu, syst=syst, weight=weights_dijuu)
                 hout['minDRuu'].fill(sample=histAxisName, channel=ch, level=lev, minDRuu=fminDRuu, syst=syst, weight=weights_dijuu)
                 hout['muu'].fill(sample=histAxisName, channel=ch, level=lev, muu=fmuu, syst=syst, weight=weights_dijuu)
                 hout['ptuu'].fill(sample=histAxisName, channel=ch, level=lev, ptuu=fptuu, syst=syst, weight=weights_dijuu)
                 hout['u0pt'].fill(sample=histAxisName, channel=ch, level=lev, u0pt=fu0pt, syst=syst, weight=weights_dijuu)
                 hout['u0eta'].fill(sample=histAxisName, channel=ch, level=lev, u0eta=fu0eta, syst=syst, weight=weights_dijuu)
   
               # Fill all the variables
               weights = weights[cut]
               fht = ht_var[cut]
               hout['counts'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut], syst=syst, weight=weights)
               hout['njets'].fill(sample=histAxisName, channel=ch, level=lev, njets=njets_var[cut], syst=syst, weight=weights)
               hout['nbtags'].fill(sample=histAxisName, channel=ch, level=lev, nbtags=nbtags_var[cut], syst=syst, weight=weights)
               hout['ht'].fill(sample=histAxisName, channel=ch, level=lev, ht=fht, syst=syst, weight=weights)
               hout['met'].fill(sample=histAxisName, channel=ch, level=lev, met=met.pt[cut], syst=syst, weight=weights)
  ### addded             
               hout['ht_atlas'].fill(sample=histAxisName, channel=ch, level=lev, ht_atlas=ht_atlas[cut], syst=syst, weight=weights)
               hout['mtw'].fill(sample=histAxisName, channel=ch, level=lev, mtw=mtw[cut], syst=syst, weight=weights) #aqui tenia flat_mtw como def en 795
               hout['met_mtw'].fill(sample=histAxisName, channel=ch, level=lev, met_mtw=met.pt[cut]+mtw[cut], syst=syst, weight=weights) #aqui tenia flat_mtw
        
               
               if lev in ['2j1b','3j1b','3j2b']:hout['DptWub'].fill(sample=histAxisName, channel=ch, level=lev, DptWub=DptWub[cut], syst=syst, weight=weights)
               if lev in ['2j1b','3j1b','3j2b']:hout['DphiWub'].fill(sample=histAxisName, channel=ch, level=lev, DphiWub=DphiWub[cut], syst=syst, weight=weights)
               if lev in ['2j1b','3j1b','3j2b']:hout['mub'].fill(sample=histAxisName, channel=ch, level=lev, mub=mub[cut], syst=syst, weight=weights)
               hout['Detalu'].fill(sample=histAxisName, channel=ch, level=lev, Detalu=Detalu[cut], syst=syst, weight=weights)
               #if lev in ['2j1b','3j1b','3j2b']:hout['Detalb'].fill(sample=histAxisName, channel=ch, level=lev, Detalb=Detalb[cut], syst=syst, weight=weights)
               #hout['DRlu'].fill(sample=histAxisName, channel=ch, level=lev, DRlu=DRlu[cut], syst=syst, weight=weights)
        
               
               if lev in ['2j1b','3j1b','3j2b']:hout['mlb'].fill(sample=histAxisName, channel=ch, level=lev, mlb=mlb[cut], syst=syst, weight=weights)
               hout['u0eta'].fill(sample=histAxisName, channel=ch, level=lev, u0eta=u0eta[cut], syst=syst, weight=weights)
               hout['absu0eta'].fill(sample=histAxisName, channel=ch, level=lev, absu0eta=abs(u0eta)[cut], syst=syst, weight=weights)              
               hout['u0pt'].fill(sample=histAxisName, channel=ch, level=lev, u0pt=u0pt[cut], syst=syst, weight=weights)
               if lev in ['2j1b','3j1b','3j2b']:hout['b0pt'].fill(sample=histAxisName, channel=ch, level=lev, b0pt=b0pt[cut], syst=syst, weight=weights)                  
               if lev in ['2j1b','3j1b','3j2b']:hout['beta'].fill(sample=histAxisName, channel=ch, level=lev, beta=b0eta[cut], syst=syst, weight=weights)
               if lev in ['2j1b','3j1b','3j2b']:hout['deltaeta'].fill(sample=histAxisName, channel=ch, level=lev, deltaeta=abs(b0eta-u0eta)[cut], syst=syst, weight=weights)
               if lev in ['2j1b','3j1b','3j2b']:hout['topmass'].fill(sample=histAxisName, channel=ch, level=lev, topmass=invariant_masses[cut], syst=syst, weight=weights)          
               
               hout['u0mass'].fill(sample=histAxisName, channel=ch, level=lev, u0mass=u0mass[cut], syst=syst, weight=weights) 
               if lev in ['2j1b','3j1b','3j2b']:hout['drub'].fill(sample=histAxisName, channel=ch, level=lev, drub=dRub[cut], syst=syst, weight=weights) 
               if lev in ['2j1b','3j1b','3j2b']:hout['absweta'].fill(sample=histAxisName, channel=ch, level=lev, absweta=w_etas[cut], syst=syst, weight=weights)
               if lev in ['2j1b','3j1b','3j2b']:hout['topeta'].fill(sample=histAxisName, channel=ch, level=lev, topeta=top_eta[cut], syst=syst, weight=weights)        
 ###


               if lev == 'incl':
                 hout['njetsnbtags'].fill(sample=histAxisName, channel=ch, level=lev, njetsnbtags=nbtagnjets[cut], syst=syst, weight=weights)
                 hout['njetsnbtags12'].fill(sample=histAxisName, channel=ch, level=lev, njetsnbtags12=nbtagnjets[cut], syst=syst, weight=weights)

               if fillAll:
                 fst = st[cut]
                 fptSumVecAll = ptSumVecAll[cut]
                 hout['st'].fill(sample=histAxisName, channel=ch, level=lev, st=fst, syst=syst, weight=weights)
                 hout['sumallpt'].fill(sample=histAxisName, channel=ch, level=lev, sumallpt=fptSumVecAll, syst=syst, weight=weights)

               if fillAll:
                 if lev != 'incl' and '0b' not in lev: # Fill jet related variables when there is at least one jet
                   fj0pt = jet0pt[cut]
                   fj0eta = jet0eta[cut]
                   jet0pt  = ak.flatten(j0.pt)
                   njetscut = njets_var[cut]
                   hout['j0pt'].fill(sample=histAxisName, channel=ch, level=lev, j0pt=fj0pt, syst=syst, weight=weights)
                   hout['j0eta'].fill(sample=histAxisName, channel=ch, level=lev, j0eta=fj0eta, syst=syst, weight=weights)
                 if lev in ['1b', '2b', '2j1b', '3j1b', '3j2b','4j1b', '4j2b', 'g5j2b']:
                   fdRlb = ak.fill_none(ak.flatten(dRlb[cut]), 0)
                   fptlb = ptSumVeclb[cut]
                   hout['ptlb'].fill(sample=histAxisName, channel=ch, level=lev, ptlb=fptlb, syst=syst, weight=weights)
                   hout['dRlb'].fill(sample=histAxisName, channel=ch, level=lev, dRlb=fdRlb, syst=syst, weight=weights)
                 if ch in ['e', 'e_fake']:
                   e = e_sel if ch == 'e' else e_fake
                   lpt  = ak.flatten(e.pt [cut])
                   leta = ak.flatten(e.eta[cut])
                   mt = ak.flatten(GetMT(e, met)[cut])
                   MVAtth = ak.flatten(e.mvaTTH[cut])
                   MVAwp = ak.flatten(e.mvaFall17V2noIso_WPL[cut])
                   hout['ept' ].fill(sample=histAxisName, channel=ch, level=lev, ept=lpt, syst=syst, weight=weights)
                   hout['eeta'].fill(sample=histAxisName, channel=ch, level=lev, eeta=leta, syst=syst, weight=weights)
                   hout['abseeta'].fill(sample=histAxisName, channel=ch, level=lev, abseeta=abs(leta), syst=syst, weight=weights)
#                   if lev in ['2j1b','3j1b','3j2b']:hout['coste'].fill(sample=histAxisName, channel=ch, level=lev, coste=cost[cut], syst=syst, weight=weights)
#                   if lev in ['2j1b','3j1b','3j2b']:hout['coste_top'].fill(sample=histAxisName, channel=ch, level=lev, coste_top=cost_top[cut], syst=syst, weight=weights)
                   
#                   if lev in ['2j1b','3j1b','3j2b']:hout['coste_1'].fill(sample=histAxisName, channel=ch, level=lev, coste_1=cos1[cut], syst=syst, weight=weights)
                   if lev in ['2j1b','3j1b','3j2b']:hout['coste_2'].fill(sample=histAxisName, channel=ch, level=lev, coste_2=cos2[cut], syst=syst, weight=weights)
#                   if lev in ['2j1b','3j1b','3j2b']:hout['coste_3'].fill(sample=histAxisName, channel=ch, level=lev, coste_3=cos2_abs[cut], syst=syst, weight=weights)
                                     
                   hout['mt'].fill(sample=histAxisName, channel=ch, level=lev, mt=mt, syst=syst, weight=weights)
                   hout['MVAtth'].fill(sample=histAxisName, channel=ch, level=lev, MVAtth=MVAtth, syst=syst, weight=weights)
                   hout['MVAwp'].fill(sample=histAxisName, channel=ch, level=lev, MVAwp=MVAwp, syst=syst, weight=weights)


                   if lev in ['2j1b','3j1b','3j2b']:               
 
 #                     mlb = (GetMlb(e[cut], goodJets[cut], int(lev[lev.find('b')-1])) )
 #                     hout['mlb'].fill(sample=histAxisName, channel=ch, level=lev, mlb=mlb, syst=syst, weight=weights)
                      if fillMVA and lev in ['2j1b','3j1b','3j2b'] and self.model_2j1b is not None and len(fht) > 0:
                        vars2j1b = [fht, fj0pt, fmjj, fdrjjmed, mlb, fdRlb, fmedDRuu, fmuu]
                        vars2j1b = [fht,abs(u0eta)[cut],fptSumVecAll,leta,fj0pt,fj0eta,u0pt[cut],ak.flatten(ptjj[cut]),fmjj,fdrjjmin,mlb[cut],ptSumVeclb[cut],
                         fdRlb,ht_atlas[cut],mtw[cut],b0eta[cut],abs(b0eta-u0eta)[cut],invariant_masses[cut],u0mass[cut],top_eta[cut],abs(w_etas)[cut],dRub[cut],cos2[cut]]
                        
                        
                        varsb10 = [abs(u0eta)[cut],abs(b0eta-u0eta)[cut],invariant_masses[cut],mt,fdrjjmed,mlb[cut],fmjj,mub[cut],fptSumVecAll,ht_atlas[cut],DptWub[cut],DphiWub[cut],Detalu[cut]]
                        


                        MVAscore_nom_pruned = self.model_pruned.predict_proba(np.column_stack(vars2j1b))[:,1] # asipara 2 clases
                        MVAscore_nom_relaxed_b10 = self.model_relaxed_b10.predict_proba(np.column_stack(varsb10))[:,1]
                     
                        hout['MVAscore_pruned'].fill(sample=histAxisName, channel=ch, level=lev, MVAscore_pruned=MVAscore_nom_pruned, syst=syst, weight=weights)
                        hout['MVAscore_relaxed_b10'].fill(sample=histAxisName, channel=ch, level=lev, MVAscore_relaxed_b10=MVAscore_nom_relaxed_b10, syst=syst, weight=weights)
     




                 elif ch in ['m', 'm_fake']:
                   m = m_sel if ch == 'm' else m_fake
                   lpt  = ak.flatten(m.pt[cut])
                   leta = ak.flatten(m.eta[cut])
                   mt = ak.flatten(GetMT(m, met)[cut])
                   hout['mpt'].fill(sample=histAxisName, channel=ch, level=lev, mpt=lpt, syst=syst, weight=weights)
                   hout['meta'].fill(sample=histAxisName, channel=ch, level=lev,meta = leta, syst=syst, weight=weights)
                   hout['absmeta'].fill(sample=histAxisName, channel=ch, level=lev,absmeta = abs(leta), syst=syst, weight=weights)
#                   if lev in ['2j1b','3j1b','3j2b']:hout['costm'].fill(sample=histAxisName, channel=ch, level=lev, costm=cost[cut], syst=syst, weight=weights)
#                   if lev in ['2j1b','3j1b','3j2b']:hout['costm_top'].fill(sample=histAxisName, channel=ch, level=lev, costm_top=cost_top[cut], syst=syst, weight=weights)
                   
#                   if lev in ['2j1b','3j1b','3j2b']:hout['costm_1'].fill(sample=histAxisName, channel=ch, level=lev, costm_1=cos1[cut], syst=syst, weight=weights)
                   if lev in ['2j1b','3j1b','3j2b']:hout['costm_2'].fill(sample=histAxisName, channel=ch, level=lev, costm_2=cos2[cut], syst=syst, weight=weights)
#                   if lev in ['2j1b','3j1b','3j2b']:hout['costm_3'].fill(sample=histAxisName, channel=ch, level=lev, costm_3=cos2_abs[cut], syst=syst, weight=weights)
                   
                   hout['mt'].fill(sample=histAxisName, channel=ch, level=lev, mt=mt, syst=syst, weight=weights)
                   if lev in ['2j1b','3j1b','3j2b']:#, '3j2b']:
 #                     mlb = (GetMlb(m[cut], goodJets[cut], int(lev[lev.find('b')-1])) )
 #                     hout['mlb'].fill(sample=histAxisName, channel=ch, level=lev, mlb=mlb, syst=syst, weight=weights)
                      # Fill MVA variables
                      if fillMVA and lev in ['2j1b','3j1b','3j2b'] and self.model_2j1b is not None and len(fht) > 0:
                        vars3j1b = [fht, fj0pt, fmjj, fdrjjmed, mlb, fdRlb, fmedDRuu, fmuu]
                        vars2j1b = [fht,abs(u0eta)[cut],fptSumVecAll,leta,fj0pt,fj0eta,u0pt[cut],ak.flatten(ptjj[cut]),fmjj,fdrjjmin,mlb[cut],ptSumVeclb[cut],
                         fdRlb,ht_atlas[cut],mtw[cut],b0eta[cut],abs(b0eta-u0eta)[cut],invariant_masses[cut],u0mass[cut],top_eta[cut],abs(w_etas)[cut],dRub[cut],cos2[cut]]
                        
                        
                        varsb10 = [abs(u0eta)[cut],abs(b0eta-u0eta)[cut],invariant_masses[cut],mt,fdrjjmed,mlb[cut],fmjj,mub[cut],fptSumVecAll,ht_atlas[cut],DptWub[cut],DphiWub[cut],Detalu[cut]]
                        
                        

                        MVAscore_nom_pruned = self.model_pruned.predict_proba(np.column_stack(vars2j1b))[:,1] # asipara 2 clases

                        MVAscore_nom_relaxed_b10 = self.model_relaxed_b10.predict_proba(np.column_stack(varsb10))[:,1]
                        
                       
                        
                        hout['MVAscore_pruned'].fill(sample=histAxisName, channel=ch, level=lev, MVAscore_pruned=MVAscore_nom_pruned, syst=syst, weight=weights)
                        hout['MVAscore_relaxed_b10'].fill(sample=histAxisName, channel=ch, level=lev, MVAscore_relaxed_b10=MVAscore_nom_relaxed_b10, syst=syst, weight=weights)                        


                 elif ch in ['em', 'ee', 'mm']:
                   llpairs = ak.combinations(l_sel[cut], 2, fields=["l0","l1"])
                   mll = (llpairs.l0+llpairs.l1).mass # Invmass for leading two leps
                   mll_flat = ak.flatten(mll)
                   lep0pt = ak.flatten(llpairs.l0.pt)
                   lep0eta = ak.flatten(llpairs.l0.eta)
                   hout['invmass'].fill(sample=histAxisName, channel=ch, level=lev, invmass=mll_flat, syst=syst, weight=weights)
                   hout['invmass2'].fill(sample=histAxisName, channel=ch, level=lev, invmass2=mll_flat, syst=syst, weight=weights)
                   hout['l0pt'].fill(sample=histAxisName, channel=ch, level=lev, lep0pt=lep0pt, syst=syst, weight=weights)
                   hout['l0eta'].fill(sample=histAxisName, channel=ch, level=lev, lep0eta=lep0eta, syst=syst, weight=weights)
                   if ch == 'ee':
                     b0 = (abs(llpairs.l0.eta) < 1.479)
                     e0 = (abs(llpairs.l0.eta) > 1.479)
                     b1 = (abs(llpairs.l1.eta) < 1.479)
                     e1 = (abs(llpairs.l1.eta) > 1.479)
                     mll_bb = mll[(b0&b1)]
                     mll_be = mll[(b0&e1)|(e0&b1)]
                     mll_ee = mll[(e0&e1)]
                     mll_bb = ak.flatten(mll_bb)
                     mll_be = ak.flatten(mll_be)
                     mll_ee = ak.flatten(mll_ee)
                     weights_bb = weights[ak.flatten(b0&b1)]
                     weights_be = weights[ak.flatten((b0&e1)|(e0&b1))]
                     weights_ee = weights[ak.flatten(e0&e1)]
                     hout['invmass_bb'].fill(sample=histAxisName, channel=ch, level=lev, invmass=mll_bb, syst=syst, weight=weights_bb)
                     hout['invmass_be'].fill(sample=histAxisName, channel=ch, level=lev, invmass=mll_be, syst=syst, weight=weights_be)
                     hout['invmass_ee'].fill(sample=histAxisName, channel=ch, level=lev, invmass=mll_ee, syst=syst, weight=weights_ee)

               # Fill scale and pdf uncertainties
               if doPDFunc and syst == 'norm':
                 scale_w = np.transpose(scaleweights[cut])*(weights)
                 pdf_w   = np.transpose(pdfweights  [cut])*(weights)
                 hout['Scales'].fill(sample=histAxisName, channel=ch, level=lev, Scales=ak.flatten(scaleweights_bins[cut]), syst="norm", weight=ak.flatten(scale_w))
                 hout['PDF']   .fill(sample=histAxisName, channel=ch, level=lev, PDF   =ak.flatten(pdfweights_bins[cut]),   syst="norm", weight=ak.flatten(pdf_w))

        return hout
  
    def postprocess(self, accumulator):
        return accumulator

