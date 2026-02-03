#!/usr/bin/env python
import lz4.frame as lz4f
import cloudpickle
import json
import pprint
import copy
import numpy as np
import awkward as ak
import coffea
import sys
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea import hist, processor
from coffea.util import load, save
from optparse import OptionParser
from coffea.analysis_tools import PackedSelection
from coffea.lumi_tools import LumiMask

from cafea.modules.GetValuesFromJsons import get_param
from cafea.analysis.objects import *
from cafea.analysis.corrections import ApplyJetSystematicsRun3, GetPUSF_run3, AttachMuonSFsRun3, AttachElecSFsRun3, AttachTrigSFsRun3, ApplyJetCorrectionsRun3
from cafea.analysis.selection import *
from cafea.modules.paths import cafea_path

'''
 Based on AN2018_210

 First:  Select a LOOSE definition for electrons and muons
 Second: Select a TIGHT definition for electrons and muons
 Third:  Choose SFs derived from TP method in Z->ll events
 
 In principle, ID variables should not depend too much on the topology of the events.
 Isolation might vary... different in DY and ttbar (phase space extrapolation bias).

 In-situ efficiency measurement ----------------------------
 Selected events --> eµ pair, OS, >=2jets, >=1btag
 Add trigger selection, overlap removal, trigger SFs, MET filters, other SFs...
 > Probe lepton: pass loose iso criteria
 > Tag   lepton: pass both loose and tight criteria
 >   Eff = Npass / Ntotal [FOR SIGNAL EVENTS]
 >   Eff = ( Npass(data) - Npass(nonprompt) ) / ( Ntot(data) - Ntot(nonprompt) )
 >   where: N(nonprompt) = (N_SS(data) - N_SS(prompt)) x R(OS/SS)
 >   Eff MC --> FROM GEN INFO + SFs from TP method in Z->ll.

Plots --> Data/MC agreement for global events
          Data/MC for passing and failing leptons, OS and SS
          Counts for efficiencies --> Eff with in-situ vs eff with "MC truth"x"Corrected"
 
'''

def isLooseElec(elec):
  return np.ones_like(elec)

def isLooseMuon(muon):
  return (muon.isGlobal)

def isMCTruthElec(elec):
  return (elec.genPartFlav == 1) | (elec.genPartFlav == 15)

def isMCTruthMuon(muon):
  return (muon.genPartFlav == 1) | (muon.genPartFlav == 15)

### GLOBAL PATHS
golden_json_path = cafea_path("data/goldenJsons/Cert_Collisions2022_356309_356615_Golden.json")

### Import yout functions to select tight and loose leptons, and overwrite the function to select MC truth if needed
ptcut = 35; eta = 2.4
isTightElec   = isElectronTight
isTightMuon   = isMuonPOGT
#isLooseElec   = 
#isLooseMuon   = 
#isMCTruthElec = 
#isMCTruthMuon = 

### Import your functions to attach the lepton SFs
AttachMuonSFs = AttachMuonSFsRun3
AttachElecSFs = AttachElecSFsRun3


class AnalysisProcessor(processor.ProcessorABC):

    global golden_json_path, isLooseElec, isLooseMuon, isTightElec, isTightMuon, isMCTruthElec, isMCTruthMuon, AttachMuonSFs, AttachElecSFs

    def __init__(self, samples):

        self._samples = samples

        # Create the histograms
        # 'name' : hist.Hist("Ytitle", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat("syst", "syst"), hist.Bin("name", "X axis (GeV)", 20, 0, 100)),
        self._accumulator = processor.dict_accumulator({

        ### Data / MC agreement, global distributions
        'counts'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('sign', 'sign'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'jet0pt'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('sign', 'sign'), hist.Bin("jet0pt",  "Leading jet $p_{T}$ (GeV)", 10, 0, 300)),
        'jet0eta'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('sign', 'sign'), hist.Bin("jet0eta", "Leading jet $\eta$ ", 12, -2.5, 2.50)),
        'jetpt'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('sign', 'sign'), hist.Bin("jetpt",  "jet $p_{T}$ (GeV)", 10, 0, 300)),
        'jeteta'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('sign', 'sign'), hist.Bin("jeteta", "jet $\eta$ ", 12, -2.5, 2.50)),
        'invmass'    : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('sign', 'sign'), hist.Bin("invmass", "$m_{\ell\ell}$ (GeV) ", 20, 0, 300)),
        'njets'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('sign', 'sign'), hist.Bin("njets",   "Jet multiplicity", 6, 0, 6)),
        "nbtags"    : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('sign', 'sign'), hist.Bin("nbtags", "Btag multiplicity ", 5, 0, 5)),
        'met'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('sign', 'sign'), hist.Bin("met",     "MET (GeV)", 10, 0, 200)),
        'ht'         : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('sign', 'sign'), hist.Bin("ht",      "H$_{T}$ (GeV)", 10, 0, 400)),

        ### Lepton distributions
        'lpt'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('lepton', 'lepton'), hist.Cat('sign', 'sign'), hist.Bin("lpt",  "Lepton $p_{T}$ (GeV)", 10, 20, 120)),
        'leta'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('lepton', 'lepton'), hist.Cat('sign', 'sign'), hist.Bin("leta", "Lepton $\eta$ ", 10, -2.5, 2.50)),
        'lpteta'       : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('lepton', 'lepton'), hist.Cat('sign', 'sign'), hist.Bin("lpt",  "Lepton $p_{T}$ (GeV)", 10, 20, 120), hist.Bin("leta", "Lepton $\eta$ ", 10, -2.5, 2.50)),

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

        ### Initialize objects
        met  = events.MET
        e    = events.Electron
        mu   = events.Muon
        jets = events.Jet

        ### Muon selection
        mu['isLoose'] = isLooseMuon(mu)
        mu['isTight'] = isTightMuon(mu)
        if not isData:
          mu['isMCtru'] = isMCTruthMuon(mu)

        ### Electron selection
        e['isLoose'] = isLooseElec(e)
        e['isTight'] = isTightElec(e)
        if not isData:
          e['isMCtru'] = isMCTruthElec(e)

        ### Attach scale factors
        if not isData:
          AttachMuonSFs(mu)
          AttachElecSFs(e)

        # Leptons
        e_sel = e[(e.isTight)]
        m_sel = mu[(mu.isTight)]
        l_sel = ak.with_name(ak.concatenate([e_sel, m_sel], axis=1), 'PtEtaPhiMCandidate')
        e_loo = e[(e.isLoose)]
        m_loo = mu[(mu.isLoose)]
        l_loo = ak.with_name(ak.concatenate([e_loo, m_loo], axis=1), 'PtEtaPhiMCandidate')

        l_sel_padded = ak.pad_none(l_sel, 2)
        l0 = l_sel_padded[:,0]
        l1 = l_sel_padded[:,1]

        leadinglep = l_sel[ak.argmax(l_sel.pt, axis=-1, keepdims=True)]
        subleadinglep = l_sel[ak.argmin(l_sel.pt, axis=-1, keepdims=True)]
        leadingpt = ak.flatten(leadinglep.pt) #ak.pad_none(l_sel.pt, 1)
        subleadingpt = ak.flatten(subleadinglep.pt) #ak.pad_none(l_sel.pt, 1)

        #### Jet selection -- cleaned with loose leptons
        jetptcut = 30
        jetptname = "pt"
        vetos_tocleanjets = ak.with_name( l_loo, "PtEtaPhiMCandidate")
        tmp = ak.cartesian([ak.local_index(jets.pt), vetos_tocleanjets.jetIdx], nested=True)
        cleanedJets = jets[~ak.any(tmp.slot0 == tmp.slot1, axis=-1)] # this line should go before *any selection*, otherwise lep.jetIdx is not aligned with the jet index
        cleanedJets["isGood"] = isTightJet(getattr(cleanedJets, jetptname), cleanedJets.eta, cleanedJets.jetId, jetPtCut=jetptcut)
        goodJets = cleanedJets[cleanedJets.isGood]

        btagwp = 0.5
        goodJets['isBtag'] = (goodJets.btagDeepFlavB > btagwp)
 
        ### Event selection ---------------------------- For control plots, based on two tight leptons
        ### Channel
        events['Eisem'] = (ak.num(m_sel) == 1) & (ak.num(e_sel) == 1)
        events['Eismm'] = (ak.num(m_sel) == 2) & (ak.num(e_sel) == 0)
        events['Eisee'] = (ak.num(m_sel) == 0) & (ak.num(e_sel) == 2)

        ### Same sign
        events['EisOS'] = (ak.prod(l_sel.charge, axis=1) == -1)
        events['EisSS'] = (ak.prod(l_sel.charge, axis=1) ==  1)

        ### Initialize weights
        if (isData): genw = np.ones_like(events["event"])
        else:        genw = events["genWeight"]
        weights_ones = np.ones_like(events["event"])
        weights_dict = coffea.analysis_tools.Weights(len(events),storeIndividual=True)
        weights_dict.add("norm",genw if isData else (xsec/sow)*genw)

        ################################ CORRECTIONS
        if not isData: # Apply SFs
          AddSFsRun3(events, l_sel)
          weights_dict.add("lepSF_muon", copy.deepcopy(events.sf_muon), copy.deepcopy(events.sf_hi_muon), copy.deepcopy(events.sf_lo_muon))
          weights_dict.add("lepSF_elec", copy.deepcopy(events.sf_elec), copy.deepcopy(events.sf_hi_elec), copy.deepcopy(events.sf_lo_elec))

        ############################### Trigger, MET filters...
        trig = trgPassNoOverlap(events,isData,dataset,year)  
        METfilters = PassMETfilters(events,isData)
        lumi_mask = np.ones_like(events['event'], dtype=bool)
        if isData and golden_json_path != '':
          lumi_mask = LumiMask(golden_json_path)(events.run,events.luminosityBlock)

        # Variables FOR CONTROL PLOTS
        counts = np.ones_like(events['event'], dtype=float)
        njets  = ak.num(goodJets)
        nbtags = ak.num(goodJets['isBtag'])
        ht = ak.sum(goodJets.pt,axis=-1)
        metpt = met.pt

        # 'j0pt', 'j0eta', 'jpt', 'jeta'
        j0 = goodJets[ak.argmax(goodJets.pt,axis=-1,keepdims=True)]

        # invmass
        llpairs = ak.combinations(l_sel, 2, fields=["l0","l1"])
        invmass = (llpairs.l0+llpairs.l1).mass
        mllvalues = np.where(ak.num(invmass)==0, [[0]], invmass)
        mllvalues = np.where(ak.num(mllvalues)>1, [[0]], mllvalues)
        mllvalues = ak.flatten(mllvalues, axis=1)

        channels = ['em', 'ee', 'mm'] 
        levels   = ['incl', 'g2jets', 'g1btag']
        signs    = ['OS', 'SS']
        baseSelect = (njets >= 0)&(trig)&(METfilters)&(lumi_mask)&(mllvalues>20)&(leadingpt>20)&(subleadingpt>20)

        selections = PackedSelection(dtype='uint64')
        selections.add("em", ( (events.Eisem)))
        selections.add("ee", ( (events.Eisee)))
        selections.add("mm", ( (events.Eismm)))
        selections.add("OS", ( (events.EisOS)))
        selections.add("SS", ( (events.EisSS)))
        selections.add("incl",  (baseSelect))
        selections.add("g1jets", (njets >= 1))
        selections.add("g2jets", (njets >= 2))
        selections.add("g1btag", (njets >= 2)&(nbtags>=1))

        ##### Loop over the hists we want to fill
        hout = self.accumulator.identity()

        ### CONTROL PLOTS
        for ch in channels:
          for lev in levels:
            for sign in signs:
              cuts = [ch] + [lev] + [sign] + ['incl']
              cut   = selections.all(*cuts)
              weights = weights_dict.weight(None)
              weight = weights[cut]

              hout['counts'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut], sign=sign, weight=weight)
              hout['njets'].fill(sample=histAxisName, channel=ch, level=lev, njets=njets[cut], sign=sign, weight=weight)
              hout['nbtags'].fill(sample=histAxisName, channel=ch, level=lev, nbtags=nbtags[cut], sign=sign, weight=weight)
              hout['ht'].fill(sample=histAxisName, channel=ch, level=lev, ht=ht[cut], sign=sign, weight=weight)
              hout['met'].fill(sample=histAxisName, channel=ch, level=lev, met=met.pt[cut], sign=sign, weight=weight)

              hout['invmass'].fill(sample=histAxisName, channel=ch, level=lev, invmass=mllvalues[cut], sign=sign, weight=weight)

              cutg1jet = selections.all(*(cuts+['g1jets']))
              jet0pt   = ak.flatten(j0.pt[cutg1jet])
              jet0eta  = ak.flatten(j0.eta[cutg1jet])
              jetpt    = ak.flatten(goodJets.pt[cutg1jet])
              jeteta   = ak.flatten(goodJets.eta[cutg1jet])
              weights_cutg1jet = weights[cutg1jet]
              hout['jet0pt'].fill(sample=histAxisName, channel=ch, level=lev, jet0pt=jet0pt, sign=sign, weight=weights_cutg1jet)
              hout['jet0eta'].fill(sample=histAxisName, channel=ch, level=lev, jet0eta=jet0eta, sign=sign, weight=weights_cutg1jet)
              
              weights_jets, _ = ak.broadcast_arrays(weights_cutg1jet, goodJets.pt[cutg1jet])
              weights_jets = ak.flatten(weights_jets)
              hout['jetpt'].fill(sample=histAxisName, channel=ch, level=lev, jetpt=jetpt, sign=sign, weight=weights_jets)
              hout['jeteta'].fill(sample=histAxisName, channel=ch, level=lev, jeteta=jeteta, sign=sign, weight=weights_jets)


        ############################
        ### TP selection ------------------------------- We select TP candidate pairs! Then treat equally all valid unique TP pairs
        llTPpairs = ak.with_name(ak.concatenate([e_loo, m_loo], axis=1), 'PtEtaPhiMCandidate')
        TPpairs1 = ak.combinations(llTPpairs, 2, fields=["tag","probe"])
        TPpairs2 = ak.combinations(llTPpairs, 2, fields=["probe","tag"])

        # Pad selection, weights, etc, to num of elements in llTPpairs...
        indices = ak.local_index(TPpairs1)
        njetsTP1 = ak.broadcast_arrays(njets, indices)[0]
        nbtagsTP1 = ak.broadcast_arrays(nbtags, indices)[0]
        njetsTP2 = ak.broadcast_arrays(njets, indices)[0]
        nbtagsTP2 = ak.broadcast_arrays(nbtags, indices)[0]

        weightsnom = genw if isData else (xsec/sow)*genw
        weightsTP1 = ak.broadcast_arrays(weightsnom, indices)[0]
        weightsTP2 = ak.broadcast_arrays(weightsnom, indices)[0]
        
        TPpairs = ak.with_name(ak.concatenate([TPpairs1, TPpairs2], axis=1), 'PtEtaPhiMCandidate')
        njetsTP = ak.with_name(ak.concatenate([njetsTP1, njetsTP2], axis=1), 'int64')
        nbtagsTP = ak.with_name(ak.concatenate([nbtagsTP1, nbtagsTP2], axis=1), 'int64')
        weightsTP = ak.with_name(ak.concatenate([weightsTP1, weightsTP2], axis=1), 'float64')

        # ... and then flatten them
        TPpairs = ak.flatten(TPpairs, axis=1)
        njetsTP = ak.flatten(njetsTP, axis=1)
        nbtagsTP = ak.flatten(nbtagsTP, axis=1)
        weightsTP = ak.flatten(weightsTP, axis=1)
        
        goodTag = (TPpairs.tag.pt>25)&(abs(TPpairs.tag.eta)<2.4)&(TPpairs.tag.isTight)
        isOS = (TPpairs.tag.charge*TPpairs.probe.charge < 0)
        isSS = (TPpairs.tag.charge*TPpairs.probe.charge > 0)
        tagElec = (np.abs(TPpairs.tag.pdgId)==11)
        tagMuon = (np.abs(TPpairs.tag.pdgId)==13)
        probeElec = (np.abs(TPpairs.probe.pdgId)==11)
        probeMuon = (np.abs(TPpairs.probe.pdgId)==13)
        isem = (tagElec&probeMuon) | (tagMuon&probeElec)
        probePass = (TPpairs.probe.isTight)
        if not isData:
          genmatch = (TPpairs.probe.isMCtru)

        base = (goodTag)&(isem)
        selectionsTP = PackedSelection(dtype='uint64')

        levels   = ['incl', 'g2jets', 'g1btag']
        signs    = ['OS', 'SS']
        channels = ['eprobe', 'muprobe'] 
        lepcat = ['pass', 'fail']

        selectionsTP.add("base", base)
        selectionsTP.add("incl", base)
        selectionsTP.add("g1btag", (njetsTP >= 2)&(nbtagsTP>=1))
        selectionsTP.add("g2jets", (njetsTP >= 2))

        selectionsTP.add("OS", isOS)
        selectionsTP.add("SS", isSS)

        selectionsTP.add("eprobe", probeElec)
        selectionsTP.add("muprobe", probeMuon)

        selectionsTP.add("pass", probePass)
        selectionsTP.add("fail", probePass==False)

        if not isData:
          lepcat += ['genmatch']
          selectionsTP.add("genmatch", genmatch)

        if not isData:
          SFTag   = TPpairs.tag.sf_nom_muon * TPpairs.tag.sf_nom_elec
          SFProbe = TPpairs.probe.sf_nom_muon * TPpairs.probe.sf_nom_elec
          SFlep = SFTag * SFProbe
          weightsTP = weightsTP * SFlep

        pt  = TPpairs.probe.pt
        eta = TPpairs.probe.eta

        for lev in levels:
          for ch in channels:
            for sign in signs:
              for lep in lepcat:
                cuts = [ch] + [lev] + [sign] + [lep] + ['base']
                cut   = selectionsTP.all(*cuts)
                weight = weightsTP[cut]

                hout['lpt'].fill(sample=histAxisName, channel=ch, level=lev, sign=sign, lepton=lep, lpt=pt[cut], weight=weight)
                hout['leta'].fill(sample=histAxisName, channel=ch, level=lev, sign=sign, lepton=lep, leta=eta[cut], weight=weight)
                hout['lpteta'].fill(sample=histAxisName, channel=ch, level=lev, sign=sign, lepton=lep, lpt=pt[cut], leta=eta[cut], weight=weight)

        return hout
  
    def postprocess(self, accumulator):
        return accumulator

