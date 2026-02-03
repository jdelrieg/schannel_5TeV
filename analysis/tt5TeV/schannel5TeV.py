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

import warnings
# ignore userwarnings and runtime warnings
warnings.filterwarnings("ignore", category=UserWarning) # evaluating the MVA because the dataset has not names for the variables...
warnings.filterwarnings("ignore", category=RuntimeWarning) # variables in the nanoAOD 

import copy

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
  cleanedJets["isGood"] = isTightJet(getattr(cleanedJets, jetptname), cleanedJets.eta, cleanedJets.jetId, jetPtCut=jetptcut)
  goodJets = cleanedJets[cleanedJets.isGood]
  goodJets["isBtag"] = (goodJets.btagDeepB > btagwp)
  btagSF = np.ones(len(goodJets), dtype=float) 
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

  if 'metg30' in cuts:
    metcat = (metpt > 30)
  elif 'metg20' in cuts:
    metcat = (metpt > 20)
  elif 'metg25' in cuts:
    metcat = (metpt > 25)
  elif 'metg15' in cuts:
    metcat = (metpt > 15)
  elif 'metl15' in cuts:
    metcat = (metpt < 15)
  elif 'metl20' in cuts:
    metcat = (metpt < 20)
  elif 'metl25' in cuts:
    metcat = (metpt < 25)

  if 'g2jets' in cuts:
    jetcat = (njets >= 2)
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
  elif '2j1b' in cuts:
    jetcat = (njets == 2) & (nbtags == 1)
  elif '3j1b' in cuts:
    jetcat = (njets == 3) & (nbtags == 1)
  elif '3j2b' in cuts:
    jetcat = (njets == 3) & (nbtags >= 2)
  elif '4j1b' in cuts:
    jetcat = (njets == 4) & (nbtags == 1)
  elif '4j2b' in cuts:
    jetcat = (njets == 4) & (nbtags >= 2)
  
  elif 'g4j1b' in cuts:
    jetcat = (njets >= 4) & (nbtags == 1)
  elif 'g4j2b' in cuts:
    jetcat = (njets >= 4) & (nbtags >= 2)

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


















class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, model):

        self._samples = samples
        self.model = model
        self.model_3j1b = self.model; self.model_3j2b = None
        if isinstance(model, list) and len(model) >= 2:
          self.model_3j1b = model[0]
          self.model_3j2b = model[1]

        # Create the histograms
        # 'name' : hist.Hist("Ytitle", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat("syst", "syst"), hist.Bin("name", "X axis (GeV)", 20, 0, 100)),
        self._accumulator = processor.dict_accumulator({
        'dummy'      : hist.Hist("Dummy", hist.Cat("sample", "sample"), hist.Bin("dummy", "Number of events", 1, 0, 1)),
        'counts'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'l0pt'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("lep0pt",  "Leading lepton $p_{T}$ (GeV)", 10, 20, 120)),
        'PDF'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("PDF",     "Counts", 33, 0, 33)),
        'Scales'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("Scales",  "Counts", 9, 0, 9)),
        'l0eta'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("lep0eta", "Leading lepton $\eta$ ", 10, -2.5, 2.50)),
        'ept'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("ept",  "Electron $p_{T}$ (GeV)", 10, 20, 120)),
        'eeta'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("eeta", "Electron $\eta$ ", 10, -2.5, 2.50)),
        'mpt'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("mpt",  "Muon $p_{T}$ (GeV)", 10, 20, 120)),
        'meta'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("meta", "Muon $\eta$ ", 10, -2.5, 2.50)),
        'j0pt'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("j0pt",  "Leading jet $p_{T}$ (GeV)", 10, 0, 300)),
        'j0eta'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("j0eta", "Leading jet $\eta$ ", 12, -2.5, 2.50)),
        'invmass'    : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("invmass", "$m_{\ell,\ell}$ (GeV) ", 20, 0, 300)),
        'invmass2'   : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("invmass2", "$m_{\ell,\ell}$ (GeV) ", 30, 70, 110)),
        'invmass_bb' : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("invmass", "$m_{\ell,\ell}$ (GeV) ", 30, 70, 110)),
        'invmass_be' : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("invmass", "$m_{\ell,\ell}$ (GeV) ", 30, 70, 110)),
        'invmass_ee' : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("invmass", "$m_{\ell,\ell}$ (GeV) ", 30, 70, 110)),
        'njets'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("njets",   "Jet multiplicity", 10, 0, 10)),
        'nbtags'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("nbtags",  "b-tag multiplicity", 4, 0, 4)),
        'met'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("met",     "MET (GeV)", 40, 0, 200)),
        'metnocut'   : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("metnocut",     "MET (GeV)", 40, 0, 200)),
        'ht'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("ht",      "H$_{T}$ (GeV)", 25, 0, 500)),
        'mt'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("mt",      "m$_{T}$ (GeV)", 10, 0, 150)),
        'mlb'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("mlb",     "m(l,b) (GeV)", 12, 0, 400)),


        'minDRjj'    : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("minDRjj", "min $\Delta$R(j,j) ", 10, 0, 3)),
        'maxDRjj'    : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("maxDRjj", "max $\Delta$R(j,j) ", 10, 1, 4)),
        'medianDRjj' : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("medianDRjj", r'$\Delta R_{med}$(j,j)', 15, 0, 4.5)),
        'mjj'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("mjj",     "m(j,j) (GeV)", 10, 0, 200)),
        'ptjj'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("ptjj",    "p$_{T}$(j,j) (GeV)", 15, 0, 300)),

        'njetsnbtags': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("njetsnbtags","(j,b)", 16, 0, 16)),
        'njetsnbtags12': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("njetsnbtags12","(j,b)", 12, 4, 16)),
  
        'u0pt'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("u0pt",  "Leading u-jet $p_{T}$ (GeV)", 10, 0, 300)),
        'u0eta'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("u0eta", "Leading u-jet $\eta$ ", 12, -2.5, 2.50)),

        'minDRuu'    : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("minDRuu", "min$\Delta$R(u,u) ", 10, 0, 3)),
        'medianDRuu' : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("medianDRuu", "median$\Delta$R(u,u) ", 15, 0, 4.5)),

        'muu'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("muu",     "m(u,u') (GeV)", 15, 10, 310)),
        'ptuu'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("ptuu",    "p$_{T}$(uu) (GeV)", 15, 0, 300)),

         # ptSumVecAll, ptSumVeclb, dRlb = GetJetLepVar(goodJets, leps)
        'st'           : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("st",      "S$_{T}$ (GeV)", 20, 0, 800)),
        'ptlb'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("ptlb",    "p$_{T}$($\ell,$b) (GeV)", 10, 0, 300)),
        'sumallpt'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("sumallpt", "$\sum_\mathrm{j,\ell}\,\mathrm{p}_{T}$ (GeV)", 10, 0, 300)),
        'dRlb'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("dRlb", "min$\Delta$R($\ell,$b) ", 10, 0, 3)),
        'MVAscore'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("MVAscore", "MVA score", 20, 0, 1)),

        'MVAtth'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("MVAtth", "MVA tth", 20, -1, 1)),
        'MVAwp'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("MVAwp", "MVA wp", 2, 0, 2)),
        'Nvtx'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("Nvtx",     "N$_{vtx}$", 12, 0, 12)),


        'counts_metg20': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_metl20': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_metg15': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_metl15': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_metg30': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_metl30': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_metg25': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'counts_metl25': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("counts",  "Counts", 1, 0, 10)),

        'pttrig' : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Bin("pttrig",  "Lepton $p_\mathrm{T}$", 20, 0, 100)),
        'etatrig': hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Bin("etatrig",  "Lepton $\eta$", [0, 1.479, 2.4])),

        # processor.column_accumulator
        # Regions: 
        #  A: 2j1b, 3j1b, 3j2b,   
        #  B: 3j1b
        # Vars:
        # njets, nbtags, ht, dRlb, st, ptlb, sumallpt, ptuu, muu, medianDRuu, u0pt, u0eta, ptjj, mjj, medianDRjj, minDRjj, mlb, mt
        '3j1b_njets' : processor.column_accumulator(np.array([])),
        '3j1b_nbtags' : processor.column_accumulator(np.array([])),
        '3j1b_ht' : processor.column_accumulator(np.array([])),
        '3j1b_st' : processor.column_accumulator(np.array([])),
        '3j1b_sumAllPt' : processor.column_accumulator(np.array([])),
        '3j1b_leta' : processor.column_accumulator(np.array([])),
        '3j1b_j0pt' : processor.column_accumulator(np.array([])),
        '3j1b_j0eta' : processor.column_accumulator(np.array([])),
        '3j1b_u0pt' : processor.column_accumulator(np.array([])),
        '3j1b_u0eta' : processor.column_accumulator(np.array([])),
        '3j1b_ptjj' : processor.column_accumulator(np.array([])),
        '3j1b_mjj' : processor.column_accumulator(np.array([])),
        '3j1b_medianDRjj' : processor.column_accumulator(np.array([])),
        '3j1b_minDRjj' : processor.column_accumulator(np.array([])),
        '3j1b_mlb' : processor.column_accumulator(np.array([])),
        '3j1b_mt' : processor.column_accumulator(np.array([])),
        '3j1b_ptsumveclb' : processor.column_accumulator(np.array([])),
        '3j1b_drlb' : processor.column_accumulator(np.array([])),
        '3j1b_druu' : processor.column_accumulator(np.array([])),
        '3j1b_druumedian' : processor.column_accumulator(np.array([])),
        '3j1b_muu' : processor.column_accumulator(np.array([])),
        '3j1b_ptuu' : processor.column_accumulator(np.array([])),

        '3j2b_njets' : processor.column_accumulator(np.array([])),
        '3j2b_nbtags' : processor.column_accumulator(np.array([])),
        '3j2b_ht' : processor.column_accumulator(np.array([])),
        '3j2b_st' : processor.column_accumulator(np.array([])),
        '3j2b_sumAllPt' : processor.column_accumulator(np.array([])),
        '3j2b_leta' : processor.column_accumulator(np.array([])),
        '3j2b_j0pt' : processor.column_accumulator(np.array([])),
        '3j2b_j0eta' : processor.column_accumulator(np.array([])),
        '3j2b_u0pt' : processor.column_accumulator(np.array([])),
        '3j2b_u0eta' : processor.column_accumulator(np.array([])),
        '3j2b_ptjj' : processor.column_accumulator(np.array([])),
        '3j2b_mjj' : processor.column_accumulator(np.array([])),
        '3j2b_medianDRjj' : processor.column_accumulator(np.array([])),
        '3j2b_minDRjj' : processor.column_accumulator(np.array([])),
        '3j2b_mlb' : processor.column_accumulator(np.array([])),
        '3j2b_mt' : processor.column_accumulator(np.array([])),
        '3j2b_ptsumveclb' : processor.column_accumulator(np.array([])),
        '3j2b_drlb' : processor.column_accumulator(np.array([])),
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
        year         = self._samples[dataset]["year"]
        xsec         = self._samples[dataset]["xsec"]
        sow          = self._samples[dataset]["nSumOfWeights"]
        isData       = self._samples[dataset]["isData"]
        isSystSample = ('mtop' in histAxisName) or ('hdamp' in histAxisName) or ('UE' in histAxisName)
        doPS         = (histAxisName in ['tt', 'ttPS']) and events.PSWeight is not None and len(events.PSWeight[0])>=4
        doPDFunc = "sumPDFWeights" in self._samples[dataset]
        '''
        PDF set: https://lhapdf.hepforge.org/pdfsets
        b = Events.GetBranch('LHEPdfWeight')
        print(b.GetTitle())
        '''
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
        PU   = events.PV.npvsGood
        
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
        #mu['isdilepton']=MuonLoose(mu.pt, mu.eta, mu.dxy, mu.dz, mu.sip3d, mu.mediumPromptId, mu.btagDeepB, ptCut=10, etaCut=2.4)

        # Electron selection
        #e['pt'] = GetElecPt(e.pt, e.eta, isdata = isData)
        #e['pt'] = GetElecPtSmear(e.pt, e.eta, isdata = isData)
        GetElecScale5TeV(e, run=306936, isData=isData)
        e['isLoose'] = ElecLoose(e.pt, e.eta, e.lostHits, e.sip3d, e.dxy, e.dz, e.btagDeepB, e.convVeto, e.mvaFall17V2noIso_WPL, 20, 2.4)
        e['isMVA']   = ElecMVA(e.miniPFRelIso_all, e.mvaTTH)
        #e['isTight'] = (e.mvaFall17V2Iso_WP80)&(e.pt>20)&(np.abs(e.eta)<2.4)
        #e['isNoIso'] = (e.mvaFall17V2noIso_WP80)&(e.pt>20)&(np.abs(e.eta)<2.4)&(e.isTight == 0)
        #e['isdilepton'] = ElecLoose(e.pt, e.eta, e.lostHits, e.sip3d, e.dxy, e.dz, e.btagDeepB, e.convVeto, e.mvaFall17V2noIso_WPL, 10, 2.5)
        

        # Build loose collections
        m_sel = mu[mu.isLoose & mu.isMVA]
        e_sel = e[e.isLoose & e.isMVA]
        m_fake = mu[mu.isLoose & (mu.isMVA == 0)]
        e_fake = e[e.isLoose & (e.isMVA == 0)]
        
        #e_dilep =e[e.isdilepton & e.isMVA]
        #m_dilep =mu[mu.isdilepton & mu.isMVA]
        #m_sel = mu[mu.isTight & mu.isIso]
        #m_fake = mu[mu.isTight & (mu.isIso == 0)]
        #e_sel = e[e.isTight]
        #e_fake = e[e.isNoIso]
        e0 = e_sel[ak.argmax(e_sel.pt,axis=-1,keepdims=True)]
        m0 = m_sel[ak.argmax(m_sel.pt,axis=-1,keepdims=True)]
  
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
        events['ise' ] = (ak.num(m_sel) == 0) & (ak.num(e_sel) == 1) #& (ak.num(m_dilep) == 0)
        events['ism' ] = (ak.num(m_sel) == 1) & (ak.num(e_sel) == 0) #& (ak.num(e_dilep) == 0)
        events['ise_fake' ] = (ak.num(m_sel) == 0) & ((ak.num(e_sel)) == 0) & (ak.num(e_fake) == 1)
        events['ism_fake' ] = (ak.num(m_sel) == 0) & ((ak.num(e_sel)) == 0) & (ak.num(m_fake) == 1)
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
          
          #print('Pts:',*cleanedJets.pt)
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
        
        # Recommendations for 94X: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
        # Medium DeepCSV for 94X wp = 0.4941 , tight = 0.8001, loose = 0.1522
        # Medium DeepJet for 94X wp = 0.3033 
        # We're using DeepCSV...
        wp = 0.4941#0.8001 #0.4941 # medium
        goodJets, metpt, btagSF = GetJESMETvar(corrected_jets, met, jetptcut, btagwp=wp, isData=isData, var='')
        goodJets_norm = goodJets
        metpt_norm = metpt

        ######### SFs, weights, systematics ##########
        # Btag SF following 1a) in https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods 
        if not isData:
          flav_cut=goodJets.hadronFlavour>=4       #mask to divide flavour of jets into light(udsg) and heavy(c,b)
          goodJets1=goodJets[flav_cut]
          goodJets2=goodJets[~flav_cut]
          btagSF_h, btagSF_hUp, btagSF_hDo = GetBtagSF5TeV(goodJets1.pt, goodJets1.eta, goodJets1.hadronFlavour, goodJets1.isBtag, True)
          btagSF_d, btagSF_dUp, btagSF_dDo = GetBtagSF5TeV(goodJets2.pt, goodJets2.eta, goodJets2.hadronFlavour, goodJets2.isBtag, True) #Now we take two btagSF divided into jet flavour
          
          #btagSF, btagSFUp, btagSFDo = GetBtagSF5TeV(goodJets.pt, goodJets.eta, goodJets.hadronFlavour, goodJets.isBtag, True)    #Originally taking just one btagSF



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
            weights_dict[ch_name].add("btagSF_h", ak.copy(btagSF_h), ak.copy(btagSF_hUp), ak.copy(btagSF_hDo))
            weights_dict[ch_name].add("btagSF_d", ak.copy(btagSF_d), ak.copy(btagSF_dUp), ak.copy(btagSF_dDo))
            #weights_dict[ch_name].add("btagSF", ak.copy(btagSF), ak.copy(btagSFUp), ak.copy(btagSFDo))
            weights_dict[ch_name].add("prefire", ak.copy(prefweight), ak.copy(prefweightUp), ak.copy(prefweightDown))
        
     
          # PS = ISR, FSR (on ttPS only)
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
        systJES = ['JER','MC', 'AbsStat', 'AbsScale', 'AbsMPF', 'Frag', 'ECAL', 'HCAL', 'Flavor', 'RelStat', 'RelPt', 'RelBal', 'RelJER', 'L3Res','MET_UnclusteredEnergy'] if splitJES else ['Total']
        systJets = [x + 'Up' for x in systJES] + [x + 'Down' for x in systJES]
        if not isData and not isSystSample: systList = systList + ["elecSFUp","elecSFDown", "muonSFUp", "muonSFDown", "btagSF_hUp", "btagSF_hDown", "btagSF_dUp", "btagSF_dDown", "prefireUp", "prefireDown", "trigSFUp", "trigSFDown"]+systJets#, "trigSFUp", "trigSFDown"] + systJets
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
        
       

        selections.add("incl", ak.ones_like(met.pt, dtype=bool))
        metfilters = PassMETfilters(events,isData)
        selections.add("METfilters", metfilters)
        # Counts
        counts = np.ones_like(events['event'], dtype=float)
 
        # Initialize the out object
        hout = self.accumulator.identity()
        channels = ['e', 'm', 'e_fake', 'm_fake'] #['em', 'e', 'm', 'ee', 'mm', 'e_fake', 'm_fake'] 
    
        #levels = ['incl', 'g2jets', 'g3jets', 'g4jets', '2j1b', '3j1b', '3j2b','4j1b', '4j2b', 'g5j1b', 'g5j2b']
        #levels = ['incl', 'g2jets', 'g3jets', 'g4jets', '0b','3j0b','4j0b', 'g5j0b', '2j1b', '3j1b', '3j2b','4j1b', '4j2b', 'g5j1b', 'g5j2b']
        levels = ['3j1b','3j2b','g4j1b','g4j2b','4j1b','4j2b','g5j1b','g5j2b','g3jets']
   

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
        j0_nom, drjj_nom, drjjmedian_nom, mjj_nom, ptjj_nom,drjjmax_nom = GetJetVariables(goodJets)
        u0_nom, druu_nom, druumedian_nom, muu_nom, ptuu_nom,druumax_nom = GetJetVariables(goodJets[(goodJets.isBtag==0)])
        ptSumVecAll_nom, ptSumVeclb_nom, dRlb_nom, st_nom   = GetJetLepVar(goodJets, leps)
        ptSumVecAll_fak, ptSumVeclb_fak, dRlb_fak, st_fak   = GetJetLepVar(goodJets, fakes)


        # Loop over the hists we want to fill
        for syst in systList:
          j0, drjj, drjjmedian, mjj, ptjj, drjjmax   = (j0_nom, drjj_nom, drjjmedian_nom, mjj_nom, ptjj_nom,drjjmax_nom)
          u0, druu, druumedian, muu, ptuu, druumax   = (u0_nom, druu_nom, druumedian_nom, muu_nom, ptuu_nom, druumax_nom)
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
            j0, drjj, drjjmedian, mjj, ptjj, drjjmax   = GetJetVariables(goodJets)
            u0, druu, druumedian, muu, ptuu, druumax   = GetJetVariables(goodJets[(goodJets.isBtag==0)])
            ptSumVecAll, ptSumVeclb, dRlb, st = GetJetLepVar(goodJets, leps)

          jet0pt  = ak.flatten(j0.pt)
          jet0eta = ak.flatten(j0.eta)
          u0pt  = ak.flatten(u0.pt)
          u0eta = ak.flatten(u0.eta)

          for ch in channels:
            if syst in ['elecSFUp', 'elecSFDown', 'muonSFUp', 'muonSFDown'] and ch in ['ee', 'mm', 'em']: continue
            if ch in ['e_fake', 'm_fake']:
              ptSumVecAll, ptSumVeclb, dRlb, st = (ptSumVecAll_fak, ptSumVeclb_fak, dRlb_fak, st_fak)
            for lev in levels:
              cutschan  = [ch] + ['incl'] + ['METfilters']
              cutsel = selections.all(*cutschan)

              #cuts      = [ch] + [lev] +(['metg30'] if lev in ['2j1b', '3j1b', '3j2b'] else [])
              cuts      = [ch] + [lev] + ['metg30'] #quitar: corte en met
              cutjets = GetCutJets(cuts, syst, met_pt, njets_var, nbtags_var)
              cut = (cutsel) & (cutjets)
              cut = np.array(cut, dtype=bool)

              cutsnoMET = [ch] + [lev]
              cutjetsnomet = GetCutJets(cutsnoMET, syst, met_pt, njets_var, nbtags_var, nujets_var)
              cutnomet = (cutsel) & (cutjetsnomet)
              cutnomet = np.array(cutnomet, dtype=bool)

              weights = weights_dict[ch if not 'fake' in ch else ch[0]].weight(syst if not syst in (['norm']+systJets) else None)
             

              if fillAcc and not isData and syst=='norm' and ch in ['e', 'm'] and lev in ['3j2b', '3j1b']:
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
               hout[f'{lev}_mlb'] = processor.column_accumulator( (GetMlb(l_sel[cut], goodJets[cut], int(lev[lev.find('b')-1]))) .to_numpy())
               hout[f'{lev}_mt'] = processor.column_accumulator( ak.flatten(GetMT(l_sel, met)[cut]) .to_numpy())
               hout[f'{lev}_ptsumveclb'] = processor.column_accumulator( ptSumVeclb[cut].to_numpy())
               hout[f'{lev}_drlb'] = processor.column_accumulator( ak.flatten(dRlb[cut]).to_numpy())
               if lev == '3j1b':
                 hout[f'{lev}_druu'] = processor.column_accumulator( ak.flatten(druu[cut]).to_numpy())
                 hout[f'{lev}_druumedian'] = processor.column_accumulator( ak.flatten(druumedian[cut]).to_numpy())
                 hout[f'{lev}_muu'] = processor.column_accumulator( ak.flatten(muu[cut]).to_numpy())
                 hout[f'{lev}_ptuu'] = processor.column_accumulator( ak.flatten(ptuu[cut]).to_numpy())

              # Fill met norm histos
              cut_metg25 = (cutnomet) & GetCutJets(['metg25'], syst, met_pt, njets_var, nbtags_var)
              cut_metl25 = (cutnomet) & GetCutJets(['metl25'], syst, met_pt, njets_var, nbtags_var)
              cut_metg20 = (cutnomet) & GetCutJets(['metg20'], syst, met_pt, njets_var, nbtags_var)
              cut_metl20 = (cutnomet) & GetCutJets(['metl20'], syst, met_pt, njets_var, nbtags_var)
              cut_metg15 = (cutnomet) & GetCutJets(['metg15'], syst, met_pt, njets_var, nbtags_var)
              cut_metl15 = (cutnomet) & GetCutJets(['metl15'], syst, met_pt, njets_var, nbtags_var)
              weights_metg20 = weights[cut_metg20]; weights_metl20 = weights[cut_metl20]
              weights_metg15 = weights[cut_metg15]; weights_metl15 = weights[cut_metl15]
              weights_metg25 = weights[cut_metg25]; weights_metl25 = weights[cut_metl25]
              hout['counts_metg20'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut_metg20], syst=syst, weight=weights_metg20)
              hout['counts_metl20'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut_metl20], syst=syst, weight=weights_metl20)
              hout['counts_metg15'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut_metg15], syst=syst, weight=weights_metg15)
              hout['counts_metl15'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut_metl15], syst=syst, weight=weights_metl15)
              hout['counts_metg25'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut_metg25], syst=syst, weight=weights_metg25)
              hout['counts_metl25'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut_metl25], syst=syst, weight=weights_metl25)

              ### met dist with no cut
              hout['metnocut'].fill(sample=histAxisName, channel=ch, level=lev, metnocut=met.pt[cutnomet], syst=syst, weight=weights[cutnomet])
  
              ### We need to have 2 jets in order to calculate dijet observables
              dijet_cut = (cut) & (GetCutJets(['g2jets'], syst, met_pt, njets_var, nbtags_var, nujets_var))
              dijuu_cut = (cut) & (GetCutJets(['g2ujets'], syst, met_pt, njets_var, nbtags_var, nujets_var))
              dijet_cut = np.array(dijet_cut, dtype=bool)
              dijuu_cut = np.array(dijuu_cut, dtype=bool)
              weights_dijet = weights[dijet_cut]
              weights_dijuu = weights[dijuu_cut]
              fdrjjmed = ak.flatten(drjjmedian[dijet_cut])
              fdrjjmax = ak.flatten(drjjmax[dijet_cut])
              fdrjjmin = ak.flatten(drjj[dijet_cut])
              hout['medianDRjj'].fill(sample=histAxisName, channel=ch, level=lev, medianDRjj=fdrjjmed, syst=syst, weight=weights_dijet)
              hout['maxDRjj'].fill(sample=histAxisName, channel=ch, level=lev, maxDRjj=fdrjjmax, syst=syst, weight=weights_dijet)
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
              hout['Nvtx'].fill(sample=histAxisName, channel=ch, level=lev, Nvtx=PU[cut], syst=syst, weight=weights)

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
                  hout['mt'].fill(sample=histAxisName, channel=ch, level=lev, mt=mt, syst=syst, weight=weights)
                  hout['MVAtth'].fill(sample=histAxisName, channel=ch, level=lev, MVAtth=MVAtth, syst=syst, weight=weights)
                  hout['MVAwp'].fill(sample=histAxisName, channel=ch, level=lev, MVAwp=MVAwp, syst=syst, weight=weights)


                  if lev in ['3j1b']:
                     mlb = (GetMlb(e[cut], goodJets[cut], int(lev[lev.find('b')-1])) )
                     hout['mlb'].fill(sample=histAxisName, channel=ch, level=lev, mlb=mlb, syst=syst, weight=weights)
                     if fillMVA and lev == '3j1b' and self.model_3j1b is not None and len(fht) > 0:
                       vars3j1b = [fht, fj0pt, fmjj, fdrjjmed, mlb, fdRlb, fmedDRuu, fmuu]
                       
                      

                       MVAscore_nom_3j1b = self.model_3j1b.predict_proba(np.column_stack(vars3j1b))[:,1]
                       hout['MVAscore'].fill(sample=histAxisName, channel=ch, level=lev, MVAscore=MVAscore_nom_3j1b, syst=syst, weight=weights)

                elif ch in ['m', 'm_fake']:
                  m = m_sel if ch == 'm' else m_fake
                  lpt  = ak.flatten(m.pt[cut])
                  leta = ak.flatten(m.eta[cut])
                  mt = ak.flatten(GetMT(m, met)[cut])
                  hout['mpt'].fill(sample=histAxisName, channel=ch, level=lev, mpt=lpt, syst=syst, weight=weights)
                  hout['meta'].fill(sample=histAxisName, channel=ch, level=lev,meta = leta, syst=syst, weight=weights)
                  hout['mt'].fill(sample=histAxisName, channel=ch, level=lev, mt=mt, syst=syst, weight=weights)
                  if lev in ['3j1b']:#, '3j2b']:
                     mlb = (GetMlb(m[cut], goodJets[cut], int(lev[lev.find('b')-1])) )
                     hout['mlb'].fill(sample=histAxisName, channel=ch, level=lev, mlb=mlb, syst=syst, weight=weights)
                     # Fill MVA variables
                     if fillMVA and lev == '3j1b' and self.model_3j1b is not None and len(fht) > 0:
                       vars3j1b = [fht, fj0pt, fmjj, fdrjjmed, mlb, fdRlb, fmedDRuu, fmuu]
                       MVAscore_nom_3j1b = self.model_3j1b.predict_proba(np.column_stack(vars3j1b))[:,1]
                       hout['MVAscore'].fill(sample=histAxisName, channel=ch, level=lev, MVAscore=MVAscore_nom_3j1b, syst=syst, weight=weights)

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

