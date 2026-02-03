import numpy as np
import awkward as ak

def isLooseElec(elec):
  return np.ones_like(elec)

def isLooseMuon(muon):
  return (muon.isGlobal)

def isMCTruthElec(elec):
  return (elec.genPartFlav == 1) | (elec.genPartFlav == 15)

def isMCTruthMuon(muon):
  return (muon.genPartFlav == 1) | (muon.genPartFlav == 15)
