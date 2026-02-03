# Module to calculate the prefire weights
# Based on: https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/common/PrefireCorr.py
from tkinter import W
import uproot
from coffea import hist, lookup_tools
import os, sys
from cafea.modules.paths import cafea_path
import numpy as np
import awkward as ak
import gzip, pickle, json
from coffea.jetmet_tools import FactorizedJetCorrector, JetCorrectionUncertainty
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory, CorrectedMETFactory
from coffea.btag_tools.btagscalefactor import BTagScaleFactor
from cafea.plotter.plotter import GetHisto, GetSFfromCountsHisto, DrawEff, DrawEff2D, GetH2DfromXY

# Load the SFs and get the lookup tables

ext_prefire = lookup_tools.extractor()
ext_prefire.add_weight_sets(["pref_jetpt L1prefiring_jetpt_2017G %s"%cafea_path('data/prefire/L1prefiring_jetpt_2017G.root')])
ext_prefire.add_weight_sets(["pref_photpt L1prefiring_photonpt_2017G %s"%cafea_path('data/prefire/L1prefiring_photonpt_2017G.root')])
ext_prefire.add_weight_sets(["pref_jetpt_er L1prefiring_jetpt_2017G_error %s"%cafea_path('data/prefire/L1prefiring_jetpt_2017G.root')])
ext_prefire.add_weight_sets(["pref_photpt_er L1prefiring_photonpt_2017G_error %s"%cafea_path('data/prefire/L1prefiring_photonpt_2017G.root')])


ext_prefire.finalize()
prefEval = ext_prefire.make_evaluator()

# Min/Max Values may need to be fixed for new maps


def GetPrefireProb(pt, eta, part='jet', var=0):
    # Limits in pt and eta for jets and photons
    JetMinPt = 20; JetMaxPt = 500; JetMinEta = 2.0; JetMaxEta = 3.0
    PhotonMinPt = 20; PhotonMaxPt = 500; PhotonMinEta = 2.0; PhotonMaxEta = 3.0
    maxpt = JetMaxPt if part=='jet' else PhotonMaxPt
    minpt = JetMinPt if part=='jet' else PhotonMinPt
    maxeta = JetMaxEta if part=='jet' else PhotonMaxEta
    mineta = JetMinEta if part=='jet' else PhotonMinEta
    # Histograms contain weights as a function of abs(eta) and pt
    eta = np.abs(eta)
    # Selection of jets/photons in the correct pt and eta range
    mask = (pt > minpt) & (pt < maxpt) & (eta > mineta) & (eta < maxeta)

    pt = np.where(pt>JetMaxPt, JetMaxPt-0.1, pt)
    prob = prefEval['pref_jetpt'    if part == 'jet' else 'pref_photpt'   ](np.abs(eta), pt)
    stat = prefEval['pref_jetpt_er' if part == 'jet' else 'pref_photpt_er'](np.abs(eta), pt)
    syst = 0.1 * prob # 10% syst
    prob_up = prob + np.sqrt(stat**2 + syst**2)
    prob_do = prob - np.sqrt(stat**2 + syst**2)
    prob_up = np.where(prob_up>1, 1, prob_up)
    prob_do = np.where(prob_do<0, 0, prob_do)
    if var==0:
        return prob
    elif var==1:
        return prob_up
    elif var==-1:
        return prob_do
    else:
        print('Wrong var value, returning nominal...')
    return prob


def GetCleanAndNot(part, part_toclean_index, part_index_name='jetIdx'):
    #vetos_toclean = ak.with_name(part_toclean,"PtEtaPhiMCandidate")

    # Indices of photons pointing to jets
    part_index = getattr(part, part_index_name)

    # part_to_clean_index are the indices of good jets... lets get the indices of photons that have indices pointing to good jets
    tmp = ak.cartesian([part_index, part_toclean_index], nested=True)
    mask = ak.any(tmp.slot0 == tmp.slot1, axis=-1)

    clean = part[~mask]
    noclean = part[mask]
    return clean, noclean

def GetPrefireWeights(jets, phot, elec, var=0):
    JetMinPt = 20; JetMaxPt = 500; JetMinEta = 2.0; JetMaxEta = 3.0
    PhotonMinPt = 20; PhotonMaxPt = 500; PhotonMinEta = 2.0; PhotonMaxEta = 3.0
    jets_mask = (jets.pt > JetMinPt) & (jets.pt < JetMaxPt) & (np.abs(jets.eta) > JetMinEta) & (np.abs(jets.eta) < JetMaxEta)
    phot_mask = (phot.pt > PhotonMinPt) & (phot.pt < PhotonMaxPt) & (np.abs(phot.eta) > PhotonMinEta) & (np.abs(phot.eta) < PhotonMaxEta)
    elec_mask = (elec.pt > PhotonMinPt) & (elec.pt < PhotonMaxPt) & (np.abs(elec.eta) > PhotonMinEta) & (np.abs(elec.eta) < PhotonMaxEta)
    j_idx_good = ak.local_index(jets.pt)[jets_mask]
    p_idx_good = ak.local_index(phot.pt)[phot_mask]
    jets = jets[jets_mask]
    phot = phot[phot_mask]
    elec = elec[elec_mask]
    jets['pt'] = np.where(jets.pt>JetMaxPt, JetMaxPt-0.1, jets.pt)
    phot['pt'] = np.where(phot.pt>PhotonMaxPt, PhotonMaxPt-0.1, phot.pt)
    elec['pt'] = np.where(elec.pt>PhotonMaxPt, PhotonMaxPt-0.1, elec.pt)

    # Clean photons with selected jets
    phot_clean, phot_noclean = GetCleanAndNot(phot, j_idx_good)

    # Clean electrons with selected jets... and clean photons? TODO
    elec_clean, elec_noclean = GetCleanAndNot(elec, j_idx_good)

    # Prefire for jets
    IsNotEmpty = lambda arr : ak.sum( ( ak.num(arr, axis=-1) ) ) > 0
    WeightsFromProb = lambda x: ak.prod(ak.ones_like(x) - x, axis=-1) if IsNotEmpty(x) else [1.]
    jet_prob        = GetPrefireProb(jets.pt,         jets.eta,         part='jet',  var=var)
    phot_injet_prob = GetPrefireProb(phot_noclean.pt, phot_noclean.eta, part='phot', var=var)
    elec_injet_prob = GetPrefireProb(elec_noclean.pt, elec_noclean.eta, part='phot', var=var)

    jetpf = WeightsFromProb(jet_prob)
    photpf_injet = WeightsFromProb(phot_injet_prob)
    elecpf_injet = WeightsFromProb(elec_injet_prob)
    EG_injet = ak.min([photpf_injet, elecpf_injet], axis=-2)
    prefjet = ak.min([jetpf, EG_injet], axis=-2)

    # Photons/electrons not in jets
    phot_prob = GetPrefireProb(phot_clean.pt, phot_clean.eta, part='phot', var=var)
    elec_prob = GetPrefireProb(elec_clean.pt, elec_clean.eta, part='phot', var=var)
    photpf = WeightsFromProb(phot_prob) # ak.prod(ak.ones_like(phot_prob) - phot_prob, axis=-1) if IsNotEmpty(phot_prob) else ak.Array([[1.]])
    elecpf = WeightsFromProb(elec_prob) # ak.prod(ak.ones_like(elec_prob) - elec_prob, axis=-1) if IsNotEmpty(elec_prob) else ak.Array([[1.]])
    prefEG = ak.min([(photpf), (elecpf)], axis=-2)

    pref = prefjet * prefEG
    return pref

# Test by loading a file using nano factory
if __name__ == "__main__":
  from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
  path = '/pool/phedex/userstorage/juanr/nanoAODv4/5TeV/nanoAODnoSkim/TT_TuneCP5_5p02TeV-powheg-pythia8/190827_224001/0000/myNanoProdMc5TeVMC_NANO_35.root'
  file = uproot.open(path)
  events = NanoEventsFactory.from_root(file, entry_stop=100).events()

  pref   = GetPrefireWeights(events.Jet, events.Photon, events.Electron, var=0)
  prefUp = GetPrefireWeights(events.Jet, events.Photon, events.Electron, var=1)
  prefDo = GetPrefireWeights(events.Jet, events.Photon, events.Electron, var=-1)
  print(pref)
  print(prefUp)
  print(prefDo)
