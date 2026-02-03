#!/usr/bin/env python
import lz4.frame as lz4f
import cloudpickle
import json
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
from cafea.analysis.corrections import SFevaluator, GetBTagSF, jet_factory, jet_factory_data, met_factory, GetBtagEff, AttachMuonSF, AttachElectronSF#, AttachPerLeptonFR, GetPUSF
from cafea.analysis.selection import *
from cafea.modules.paths import cafea_path

def GetElecPt(pt, eta, ecorr = 1, isdata = False):
  eta_sep = (abs(eta) < 1.479)
  fact_leta = (1.016-0.0035) if isdata else 1.005
  fact_heta = (1.052 -0.036) if isdata else 0.992
  facts_leta = np.ones_like(pt)*fact_leta
  facts_heta = np.ones_like(pt)*fact_heta
  facts = np.where(eta_sep, facts_leta, facts_heta)
  return pt*facts

def GetElecPtSmear(pt, eta, isdata = False):
  if(isdata): return pt
  mass = 91.1876
  sigma_lowEta = 1.786/mass
  sigma_highEta = 3.451/mass
  rnd_lowEta = np.random.normal(1, sigma_lowEta, len(pt))
  rnd_highEta = np.random.normal(1, sigma_lowEta, len(pt))
  eta_sep = (abs(eta) < 1.479)
  smear = np.where(eta_sep, rnd_lowEta, rnd_highEta)
  return pt*smear

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples):

        self._samples = samples

        # Create the histograms
        # 'name' : hist.Hist("Ytitle", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat("syst", "syst"), hist.Bin("name", "X axis (GeV)", 20, 0, 100)),
        self._accumulator = processor.dict_accumulator({
        'pt'  : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("val", "val"), hist.Cat("channel", "channel"), hist.Bin("pt",  "p_{T} (GeV) ", [20, 40, 120])),
        'eta' : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("val", "val"), hist.Cat("channel", "channel"), hist.Bin("eta", "eta", [0, 1.479, 2.4])),
        'pteta' : hist.Hist("Events",hist.Cat("sample", "sample"), hist.Cat("val", "val"), hist.Cat("channel", "channel"), hist.Bin("pt",  "p_{T} (GeV) ", [20, 40, 120]), hist.Bin("abseta", "eta", [0, 1.479, 2.4])),
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

        # Initialize objects
        met  = events.MET
        e    = events.Electron
        mu   = events.Muon
        tau  = events.Tau
        jets = events.Jet

        # Pre-selection (must be updated with 5TeV definitions)
        e["idEmu"] = ttH_idEmu_cuts_E3(e.hoe, e.eta, e.deltaEtaSC, e.eInvMinusPInv, e.sieie)
        e["conept"] = coneptElec(e.pt, e.mvaTTH, e.jetRelIso)
        mu["conept"] = coneptMuon(mu.pt, mu.mvaTTH, mu.jetRelIso, mu.mediumId)
        e["btagDeepB"] = ak.fill_none(e.matched_jet.btagDeepB, -99)
        mu["btagDeepB"] = ak.fill_none(mu.matched_jet.btagDeepB, -99)

        # Muon selection
        mu["isLoose"] = MuonLoose(mu.pt, mu.eta, mu.dxy, mu.dz, mu.sip3d, mu.mediumPromptId, mu.btagDeepB, ptCut=20., etaCut=2.4)
        mu["isMVA"]= MuonMVA(mu.miniPFRelIso_all, mu.mvaTTH)

        # Electron selection
        #e_pt, e_m = GetCorrectedElec(e)
        #e['pt_corr'] = e_pt
        #e['m_corr']  = e_m
        e['pt'] = GetElecPt(e.pt, e.eta, isdata = isData)
        e['pt'] = GetElecPtSmear(e.pt, e.eta, isdata = isData)
        e['isLoose'] = ElecLoose(e.pt, e.eta, e.lostHits, e.sip3d, e.dxy, e.dz, e.btagDeepB, e.convVeto, e.mvaFall17V2noIso_WPL, 20., 2.4)
        e['isMVA']   = ElecMVA(e.miniPFRelIso_all, e.mvaTTH)

        # Build loose collections
        m_sel = mu[mu.isLoose & mu.isMVA]
        e_sel = e[e.isLoose & e.isMVA]

        l_sel = ak.with_name(ak.concatenate([e_sel, m_sel], axis=1), 'PtEtaPhiMCandidate')

        events['isem'] = (ak.num(m_sel) == 1) & (ak.num(e_sel) == 1)
        events['ismm'] = (ak.num(m_sel) == 2) & (ak.num(e_sel) == 0)
        events['isee'] = (ak.num(m_sel) == 0) & (ak.num(e_sel) == 2)
        events['ise' ] = (ak.num(m_sel) == 0) & (ak.num(e_sel) == 1)
        events['ism' ] = (ak.num(m_sel) == 1) & (ak.num(e_sel) == 0)

        ### Trigger
        trigem = (events.HLT.HIMu17) | (events.HLT.HIEle15_WPLoose_Gsf)
        trigee = (events.HLT.HIEle15_WPLoose_Gsf) | (events.HLT.HIEle17_WPLoose_Gsf)
        trigmm = (events.HLT.HIMu17)
        trige  = events.HLT.HIEle20_WPLoose_Gsf
        trigm  = events.HLT.HIMu17#HIL3Mu20

        weights = np.ones_like(events["event"])
        hout = self.accumulator.identity()

        # Events containing electron and muon
        # Muon eff
        if histAxisName!='SingleMuon': 
          mden = m_sel[(events.isem)&(trige)]
          mnum = m_sel[(events.isem)&(trige)&(trigm)]
          hout['pt']   .fill(sample=histAxisName, channel='m', val='num', pt =ak.flatten(mnum.pt ))
          hout['pt']   .fill(sample=histAxisName, channel='m', val='den', pt =ak.flatten(mden.pt ))
          hout['eta']  .fill(sample=histAxisName, channel='m', val='num', eta=ak.flatten(np.abs(mnum.eta)))
          hout['eta']  .fill(sample=histAxisName, channel='m', val='den', eta=ak.flatten(np.abs(mden.eta)))
          hout['pteta'].fill(sample=histAxisName, channel='m', val='num', pt =ak.flatten(mnum.pt), abseta=ak.flatten(np.abs(mnum.eta)))
          hout['pteta'].fill(sample=histAxisName, channel='m', val='den', pt =ak.flatten(mden.pt), abseta=ak.flatten(np.abs(mden.eta)))

        # Elec eff
        if histAxisName!='HighEGJet': 
          eden = e_sel[(events.isem)&(trigm)]
          enum = e_sel[(events.isem)&(trige)&(trigm)]
          hout['pt']   .fill(sample=histAxisName, channel='e', val='num', pt =ak.flatten(enum.pt ))
          hout['pt']   .fill(sample=histAxisName, channel='e', val='den', pt =ak.flatten(eden.pt ))
          hout['eta']  .fill(sample=histAxisName, channel='e', val='num', eta=ak.flatten(np.abs(enum.eta)))
          hout['eta']  .fill(sample=histAxisName, channel='e', val='den', eta=ak.flatten(np.abs(eden.eta)))
          hout['pteta'].fill(sample=histAxisName, channel='e', val='num', pt =ak.flatten(enum.pt), abseta=ak.flatten(np.abs(enum.eta)))
          hout['pteta'].fill(sample=histAxisName, channel='e', val='den', pt =ak.flatten(eden.pt), abseta=ak.flatten(np.abs(eden.eta)))
          
        return hout

    def postprocess(self, accumulator):
        return accumulator

if __name__ == '__main__':
    # Load the .coffea files
    outpath= './coffeaFiles/'
    samples     = load(outpath+'samples.coffea')
    topprocessor = AnalysisProcessor(samples)

