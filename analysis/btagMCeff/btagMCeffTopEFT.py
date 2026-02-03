#!/usr/bin/env python
import lz4.frame as lz4f
import cloudpickle
import json
import pprint
import coffea
import numpy as np
import awkward as ak
np.seterr(divide='ignore', invalid='ignore', over='ignore')
#from coffea.arrays import Initialize # Not used and gives error
from coffea import hist, processor
from coffea.util import load, save
from optparse import OptionParser
from coffea.analysis_tools import PackedSelection

from cafea.analysis.objects import *
from cafea.analysis.selection import *

#coffea.deprecations_as_errors = True

# In the future these names will be read from the nanoAOD files

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples):
        self._samples = samples

        # Create the histograms
        # In general, histograms depend on 'sample', 'channel' (final state) and 'cut' (level of selection)
        self._accumulator = processor.dict_accumulator({
        'jetpt'  : hist.Hist("Events", hist.Cat("WP", "WP"), hist.Cat("Flav", "Flav"), hist.Bin("pt",  "Jet p_{T} (GeV) ", 40, 0, 800)),
        'jeteta' : hist.Hist("Events", hist.Cat("WP", "WP"), hist.Cat("Flav", "Flav"), hist.Bin("eta", "Jet eta", 25, -2.5, 2.5)),
        'jetpteta' : hist.Hist("Events", hist.Cat("WP", "WP"), hist.Cat("Flav", "Flav"), hist.Bin("pt",  "Jet p_{T} (GeV) ", [20, 30, 60, 120]), hist.Bin("abseta", "Jet eta", [0, 1, 1.8, 2.4])),
        'jetptetaflav' : hist.Hist("Events", hist.Cat("WP", "WP"), hist.Bin("pt",  "Jet p_{T} (GeV) ", [20, 30, 60, 120]), hist.Bin("abseta", "Jet eta", [0, 1, 1.8, 2.4]), hist.Bin("flav", "Flavor", [0, 4, 5]) ),
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
        dataset = events.metadata['dataset']
        year   = self._samples[dataset]['year']
        xsec   = self._samples[dataset]['xsec']
        sow    = self._samples[dataset]['nSumOfWeights' ]
        isData = self._samples[dataset]['isData']

        # Initialize objects
        met = events.MET
        e   = events.Electron
        mu  = events.Muon
        j   = events.Jet
        e["conept"] = coneptElec(e.pt, e.mvaTTH, e.jetRelIso)
        mu["conept"] = coneptMuon(mu.pt, mu.mvaTTH, mu.jetRelIso, mu.mediumId)
        e["btagDeepFlavB"] = ak.fill_none(e.matched_jet.btagDeepFlavB, -99)
        mu["btagDeepFlavB"] = ak.fill_none(mu.matched_jet.btagDeepFlavB, -99)
        e["idEmu"] = ttH_idEmu_cuts_E3(e.hoe, e.eta, e.deltaEtaSC, e.eInvMinusPInv, e.sieie)

        # Muon selection
        mu["isPres"] = isPresMuon(mu.dxy, mu.dz, mu.sip3d, mu.eta, mu.pt, mu.miniPFRelIso_all)
        mu["isLooseM"] = isLooseMuon(mu.miniPFRelIso_all,mu.sip3d,mu.looseId)
        mu["isFO"] = isFOMuon(mu.pt, mu.conept, mu.btagDeepFlavB, mu.mvaTTH, mu.jetRelIso, year)
        mu["isTightLep"]= tightSelMuon(mu.isFO, mu.mediumId, mu.mvaTTH)

        leading_mu = mu[ak.argmax(mu.pt,axis=-1,keepdims=True)]
        leading_mu = leading_mu[leading_mu.isTightLep]
        
        mu = mu[mu.isTightLep]
        mu_pres = mu[mu.isPres]

        # Electron selection
        e["isPres"] = isPresElec(e.pt, e.eta, e.dxy, e.dz, e.miniPFRelIso_all, e.sip3d, getattr(e,"mvaFall17V2noIso_WPL"))
        e["isLooseE"] = isLooseElec(e.miniPFRelIso_all,e.sip3d,e.lostHits)
        e["isFO"] = isFOElec(e.conept, e.btagDeepFlavB, e.idEmu, e.convVeto, e.lostHits, e.mvaTTH, e.jetRelIso, e.mvaFall17V2noIso_WP80, year)
        e["isTightLep"] = tightSelElec(e.isFO, e.mvaTTH)


        leading_e = e[ak.argmax(e.pt,axis=-1,keepdims=True)]
        leading_e = leading_e[leading_e.isTightLep]

        e  =  e[e .isTightLep]

        nElec = ak.num(e)
        nMuon = ak.num(mu)

        twoLeps   = (nElec+nMuon) == 2
        threeLeps = (nElec+nMuon) == 3
        twoElec   = (nElec == 2)
        twoMuon   = (nMuon == 2)
        e0 = e[ak.argmax(e.pt,axis=-1,keepdims=True)]
        m0 = mu[ak.argmax(mu.pt,axis=-1,keepdims=True)]

        m_fo = mu[mu.isPres & mu.isLooseM & mu.isFO]
        e_fo = e[e.isPres & e.isLooseE & e.isFO]

        m_fo['convVeto'] = ak.ones_like(m_fo.charge); 
        m_fo['lostHits'] = ak.zeros_like(m_fo.charge); 
        l_fo = ak.with_name(ak.concatenate([e_fo, m_fo], axis=1), 'PtEtaPhiMCandidate')
        l_fo_conept_sorted = l_fo[ak.argsort(l_fo.conept, axis=-1,ascending=False)]


        # Jet selection

        jetptname = 'pt_nom' if hasattr(j, 'pt_nom') else 'pt'
        j['isGood']  = isTightJet(getattr(j, jetptname), j.eta, j.jetId)

        vetos_tocleanjets = ak.with_name( l_fo, "PtEtaPhiMCandidate")
        tmp = ak.cartesian([ak.local_index(j.pt), vetos_tocleanjets.jetIdx], nested=True)
        j['isClean'] = (~ak.any(tmp.slot0 == tmp.slot1, axis=-1))#isClean(j, e, drmin=0.4)& isClean(j, mu, drmin=0.4)
        goodJets = j[(j.isClean)&(j.isGood)]
        njets = ak.num(goodJets)
        ht = ak.sum(goodJets.pt,axis=-1)
        j0 = goodJets[ak.argmax(goodJets.pt,axis=-1,keepdims=True)]

        ### We need weights for: normalization, lepSF, triggerSF, pileup, btagSF...
        weights = coffea.analysis_tools.Weights(len(events))
        weights.add('norm', np.ones_like(met.pt))

        hout = self.accumulator.identity()
        normweights = weights.weight().flatten() 

        flavSelection = {'b': (np.abs(goodJets.hadronFlavour) == 5), 'c': (np.abs(goodJets.hadronFlavour) == 4), 'l': (np.abs(goodJets.hadronFlavour) <= 3) }

        WP = {'all' : -999., 'loose': 0.0490, 'medium': 0.2783, 'tight': 0.7100}
        btagSelection = {}
        for wp, wpvals in WP.items():
          btagSelection[wp] = (goodJets.btagDeepFlavB>wpvals) 

        for jetype in ['b', 'c', 'l']:
          for wp in WP.keys():
            mask = (flavSelection[jetype])&(btagSelection[wp])
            selectjets = goodJets[mask]
            pts     = ak.flatten(selectjets.pt)
            etas    = ak.flatten(selectjets.eta)
            absetas = ak.flatten(np.abs(selectjets.eta))
            flavarray = np.zeros_like(pts) if jetype == 'l' else (np.ones_like(pts)*(4 if jetype=='c' else 5))
            weights =  np.ones_like(pts)
            hout['jetpt'].fill(WP=wp, Flav=jetype,  pt=pts, weight=weights)
            hout['jeteta'].fill(WP=wp, Flav=jetype,  eta=etas, weight=weights)
            hout['jetpteta'].fill(WP=wp, Flav=jetype,  pt=pts, abseta=absetas, weight=weights)
            hout['jetptetaflav'].fill(WP=wp, pt=pts, abseta=absetas, flav=flavarray, weight=weights)
    
        return hout


    def postprocess(self, accumulator):
        return accumulator

if __name__ == '__main__':
    # Load the .coffea files
    outpath= './coffeaFiles/'
    samples     = load(outpath+'samples.coffea')
    topprocessor = AnalysisProcessor(samples)

