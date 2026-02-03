'''
 selection.py

 This script contains several functions that implement the some event selection. 
 The functinos defined here can be used to define a selection, signal/control region, etc.
 The functions are called with (jagged)arrays as imputs plus some custom paramenters and return a boolean mask.

'''

import numpy as np
import awkward as ak


# The datasets we are using, and the triggers in them
dataset_dict = {

    "2016" : {
        "SingleMuon" : [
            "IsoMu24",
            "IsoMu27",
        ],
        "SingleElectron" : [
            'Ele27_WPTight_Gsf'
        ],
        "DoubleMuon" : [
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
            "TripleMu_12_10_5",
        ],
        "DoubleEG" : [
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
        ],
        "MuonEG" : [
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu8_DiEle12_CaloIdL_TrackIdL",
            "DiMu9_Ele9_CaloIdL_TrackIdL",
        ]
    },

    "2017" : {
        "SingleMuon" : [
            "IsoMu24",
            "IsoMu27",
        ],
        "SingleElectron" : [
            "Ele32_WPTight_Gsf",
            "Ele35_WPTight_Gsf",
        ],
        "DoubleMuon" : [
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            "TripleMu_12_10_5",
        ],
        "DoubleEG" : [
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
        ],
        "MuonEG" : [
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu8_DiEle12_CaloIdL_TrackIdL",
            "Mu8_DiEle12_CaloIdL_TrackIdL_DZ", # Note: Listed in Andrew's thesis, but not TOP-19-001 AN
            "DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
        ]
    },

    "2018" : {
        "SingleMuon" : [
            "IsoMu24",
            "IsoMu27",
        ],
        "EGamma" : [
            "Ele32_WPTight_Gsf",
            "Ele35_WPTight_Gsf",
        ],
        "DoubleMuon" : [
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            "TripleMu_12_10_5",
        ],
        "DoubleEG" : [
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
        ],
        "MuonEG" : [
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu8_DiEle12_CaloIdL_TrackIdL",
            "Mu8_DiEle12_CaloIdL_TrackIdL_DZ",
            "DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
        ]
    },

    "2022" : {
        "SingleMuon" : [
            "IsoMu24",
        ],
        "DoubleMuon" : [
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8" 
        ],
        "EGamma" : [
            'Ele32_WPTight_Gsf',
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
            "DoubleEle25_CaloIdL_MW",
        ],
        "MuonEG" : [
            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
        ],
        "Muon" : [
            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            "IsoMu24",
        ]
    },
}


# Hard coded dictionary for figuring out overlap...
#   - No unique way to do this
#   - Note: In order for this to work properly, you should be processing all of the datastes to be used in the analysis
#   - Otherwise, you may be removing events that show up in other datasets you're not using
exclude_dict = {
    "2016": {
        "DoubleMuon"     : [],
        "DoubleEG"       : dataset_dict["2016"]["DoubleMuon"],
        "MuonEG"         : dataset_dict["2016"]["DoubleMuon"] + dataset_dict["2016"]["DoubleEG"],
        "SingleMuon"     : dataset_dict["2016"]["DoubleMuon"] + dataset_dict["2016"]["DoubleEG"] + dataset_dict["2016"]["MuonEG"],
        "SingleElectron" : dataset_dict["2016"]["DoubleMuon"] + dataset_dict["2016"]["DoubleEG"] + dataset_dict["2016"]["MuonEG"] + dataset_dict["2016"]["SingleMuon"],
    },
    "2017": {
        "DoubleMuon"     : [],
        "DoubleEG"       : dataset_dict["2017"]["DoubleMuon"],
        "MuonEG"         : dataset_dict["2017"]["DoubleMuon"] + dataset_dict["2017"]["DoubleEG"],
        "SingleMuon"     : dataset_dict["2017"]["DoubleMuon"] + dataset_dict["2017"]["DoubleEG"] + dataset_dict["2017"]["MuonEG"],
        "SingleElectron" : dataset_dict["2017"]["DoubleMuon"] + dataset_dict["2017"]["DoubleEG"] + dataset_dict["2017"]["MuonEG"] + dataset_dict["2017"]["SingleMuon"],
    },
    "2018": {
        "DoubleMuon"     : [],
        "EGamma"         : dataset_dict["2018"]["DoubleMuon"],
        "MuonEG"         : dataset_dict["2018"]["DoubleMuon"] + dataset_dict["2018"]["EGamma"],
        "SingleMuon"     : dataset_dict["2018"]["DoubleMuon"] + dataset_dict["2018"]["EGamma"] + dataset_dict["2018"]["MuonEG"],
    },
    #"2022": {
    #    "Muon"           : [],
    #    "SingleMuon"     : [],
    #    "DoubleMuon"     : dataset_dict["2022"]["SingleMuon"],
    #    "EGamma"         : dataset_dict["2022"]["SingleMuon"]+dataset_dict["2022"]["Muon"]+dataset_dict["2022"]["DoubleMuon"],
    #    "MuonEG"         : dataset_dict["2022"]["SingleMuon"]+dataset_dict["2022"]["Muon"]+dataset_dict["2022"]["DoubleMuon"] + dataset_dict["2022"]["EGamma"]
        #"SingleElectron" : dataset_dict["2016"]["DoubleMuon"] + dataset_dict["2016"]["DoubleEG"] + dataset_dict["2016"]["MuonEG"] + dataset_dict["2016"]["SingleMuon"],},
    "2022": {
        "MuonEG"         : [],
        "DoubleMuon"     : dataset_dict["2022"]["MuonEG"],
        "Muon"           : dataset_dict["2022"]["MuonEG"],
        "EGamma"         : dataset_dict["2022"]["MuonEG"] + dataset_dict["2022"]["DoubleMuon"] + dataset_dict["2022"]["Muon"],
        "SingleMuon"     : dataset_dict["2022"]["MuonEG"] + dataset_dict["2022"]["DoubleMuon"] + dataset_dict["2022"]["EGamma"]
    }
}

trigttbar = {
  "2018" : {
    "m"  : ['IsoMu24'],
    "e"  : ['Ele32_WPTight_Gsf'],
    "em" : ['Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ',
            'Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL'],
    "ee" : ['Ele23_Ele12_CaloIdL_TrackIdL_IsoVL',
            'DoubleEle25_CaloIdL_MW'],
            #'Ele32_WPTight_Gsf'],
    "mm" : ['Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8'],
  }
}

# This is a helper function called by trgPassNoOverlap
#   - Takes events objects, and a lits of triggers
#   - Returns an array the same length as events, elements are true if the event passed at least one of the triggers and false otherwise
def passesTrgInLst(events,trg_name_lst):
    tpass = np.zeros_like(np.array(events.MET.pt), dtype=np.bool)
    trg_info_dict = events.HLT

    # "fields" should be list of all triggers in the dataset
    common_triggers = set(trg_info_dict.fields) & set(trg_name_lst)

    # Check to make sure that at least one of our specified triggers is present in the dataset
    if len(common_triggers) == 0 and len(trg_name_lst):
        raise Exception("No triggers from the sample matched to the ones used in the analysis.")

    for trg_name in common_triggers:
        tpass = tpass | trg_info_dict[trg_name]
    return tpass

# This is what we call from the processor
#   - Returns an array the len of events
#   - Elements are false if they do not pass any of the triggers defined in dataset_dict
#   - In the case of data, events are also false if they overlap with another dataset
def trgPassNoOverlap(events,is_data,dataset,year):
    
    # The trigger for 2016 and 2016APV are the same
    if year == "2016APV":
        year = "2016"

    # Initialize ararys and lists, get trg pass info from events
    trg_passes    = np.zeros_like(np.array(events.MET.pt), dtype=np.bool) # Array of False the len of events
    trg_overlaps  = np.zeros_like(np.array(events.MET.pt), dtype=np.bool) # Array of False the len of events
    trg_info_dict = events.HLT
    full_trg_lst  = []

    # Get the full list of triggers in all datasets
    for dataset_name in dataset_dict[year].keys():
        full_trg_lst = full_trg_lst + dataset_dict[year][dataset_name]

    # Check if events pass any of the triggers
    trg_passes = passesTrgInLst(events,full_trg_lst)

    # In case of data, check if events overlap with other datasets
    if is_data:
        trg_passes = passesTrgInLst(events,dataset_dict[year][dataset])
        trg_overlaps = passesTrgInLst(events, exclude_dict[year][dataset])

    # Return true if passes trg and does not overlap
    return (trg_passes & ~trg_overlaps)

def PassMETfilters(events,isData):
    filter_flags = events.Flag
    filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.HBHENoiseFilter & filter_flags.HBHENoiseIsoFilter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & (isData | filter_flags.eeBadScFilter)
    return filters

# Add SFs for tt 5 TeV
def AddSFs(events, leps):
  padded_leps_1 = ak.pad_none(leps, 1)
  padded_leps_2 = ak.pad_none(leps, 2)
  events['sf_2l']    = ak.fill_none(padded_leps_2[:,0].sf_nom*padded_leps_2[:,1].sf_nom, 1)
  events['sf_2l_hi'] = ak.fill_none(padded_leps_2[:,0].sf_hi*padded_leps_2[:,1].sf_hi, 1)
  events['sf_2l_lo'] = ak.fill_none(padded_leps_2[:,0].sf_lo*padded_leps_2[:,1].sf_lo, 1)
  events['sf_1l']    = ak.fill_none(padded_leps_1[:,0].sf_nom, 1)
  events['sf_1l_hi'] = ak.fill_none(padded_leps_1[:,0].sf_hi, 1)
  events['sf_1l_lo'] = ak.fill_none(padded_leps_1[:,0].sf_lo, 1)

def PadSFs2leps(events, leps, name):
  padded_leps_2 = ak.pad_none(leps, 2)
  events[name      ] = ak.fill_none(padded_leps_2[:,0].sf_nom, 1)*ak.fill_none(padded_leps_2[:,1].sf_nom, 1)
  events[name+'_hi'] = ak.fill_none(padded_leps_2[:,0].sf_hi, 1) *ak.fill_none(padded_leps_2[:,1].sf_hi, 1)
  events[name+'_lo'] = ak.fill_none(padded_leps_2[:,0].sf_lo, 1) *ak.fill_none(padded_leps_2[:,1].sf_lo, 1)

def AddSFsRun3(events, leps):
  padded_FOs = ak.pad_none(leps, 2)
  events['sf_muon'] = ak.fill_none(padded_FOs[:,0].sf_nom_muon*padded_FOs[:,1].sf_nom_muon, 1)
  events['sf_elec'] = ak.fill_none(padded_FOs[:,0].sf_nom_elec*padded_FOs[:,1].sf_nom_elec, 1)
  events['sf_hi_muon'] = ak.fill_none(padded_FOs[:,0].sf_hi_muon*padded_FOs[:,1].sf_hi_muon, 1)
  events['sf_hi_elec'] = ak.fill_none(padded_FOs[:,0].sf_hi_elec*padded_FOs[:,1].sf_hi_elec, 1)
  events['sf_lo_muon'] = ak.fill_none(padded_FOs[:,0].sf_lo_muon*padded_FOs[:,1].sf_lo_muon, 1)
  events['sf_lo_elec'] = ak.fill_none(padded_FOs[:,0].sf_lo_elec*padded_FOs[:,1].sf_lo_elec, 1)

# 2l selection (we do not make the ss requirement here)
def add2lMaskAndSFs(events, year, isData, sampleType):
    # FOs and padded FOs
    FOs = events.l_fo_conept_sorted
    padded_FOs = ak.pad_none(FOs,2)

    # Filters and cleanups
    filter_flags = events.Flag
    filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.HBHENoiseFilter & filter_flags.HBHENoiseIsoFilter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & ((year == "2016") | filter_flags.ecalBadCalibFilter) & (isData | filter_flags.eeBadScFilter)
    cleanup = events.minMllAFAS > 12
    muTightCharge = ((abs(padded_FOs[:,0].pdgId)!=13) | (padded_FOs[:,0].tightCharge>=1)) & ((abs(padded_FOs[:,1].pdgId)!=13) | (padded_FOs[:,1].tightCharge>=1))

    # Zee veto
    Zee_veto = (abs(padded_FOs[:,0].pdgId) != 11) | (abs(padded_FOs[:,1].pdgId) != 11) | ( abs ( (padded_FOs[:,0]+padded_FOs[:,1]).mass -91.2) > 10)

    # IDs
    eleID1 = (abs(padded_FOs[:,0].pdgId)!=11) | ((padded_FOs[:,0].convVeto != 0) & (padded_FOs[:,0].lostHits==0) & (padded_FOs[:,0].tightCharge>=2))
    eleID2 = (abs(padded_FOs[:,1].pdgId)!=11) | ((padded_FOs[:,1].convVeto != 0) & (padded_FOs[:,1].lostHits==0) & (padded_FOs[:,1].tightCharge>=2))

    # 2l requirements:
    exclusive = ak.num( FOs[FOs.isTightLep],axis=-1)<3
    dilep = (ak.num(FOs)) >= 2
    pt2515 = (ak.any(FOs[:,0:1].conept > 25.0, axis=1) & ak.any(FOs[:,1:2].conept > 15.0, axis=1))
    mask = (filters & cleanup & dilep & pt2515 & exclusive & Zee_veto & eleID1 & eleID2 & muTightCharge)
    
    # MC matching requirement (already passed for data)
    if sampleType == 'prompt':
        lep1_match=((padded_FOs[:,0].genPartFlav==1) | (padded_FOs[:,0].genPartFlav == 15))    
        lep2_match=((padded_FOs[:,1].genPartFlav==1) | (padded_FOs[:,1].genPartFlav == 15))
        mask = mask & lep1_match & lep2_match
    elif sampleType =='conversions':
        lep1_match=(padded_FOs[:,0].genPartFlav==22)
        lep2_match=(padded_FOs[:,1].genPartFlav==22)
        mask = mask & ( lep1_match | lep2_match ) 
    elif sampleType == 'nonprompt':
        lep1_match=((padded_FOs[:,0].genPartFlav!=1) & (padded_FOs[:,0].genPartFlav != 15) & (padded_FOs[:,0].genPartFlav != 22))
        lep2_match=((padded_FOs[:,1].genPartFlav!=1) & (padded_FOs[:,1].genPartFlav != 15) & (padded_FOs[:,1].genPartFlav != 22))
        mask = mask & ( lep1_match | lep2_match ) 
    elif sampleType == "data":
        pass
    else:
        raise Exception(f"Error: Unknown sampleType {sampleType}.")

    events['is2l'] = ak.fill_none(mask,False)

    # SFs
    events['sf_2l'] = padded_FOs[:,0].sf_nom*padded_FOs[:,1].sf_nom
    events['sf_2l_hi'] = padded_FOs[:,0].sf_hi*padded_FOs[:,1].sf_hi
    events['sf_2l_lo'] = padded_FOs[:,0].sf_lo*padded_FOs[:,1].sf_lo

    # SR:
    events['is2l_SR'] = (padded_FOs[:,0].isTightLep) & (padded_FOs[:,1].isTightLep)
    events['is2l_SR'] = ak.fill_none(events['is2l_SR'],False)



# 3l selection
def add3lMaskAndSFs(events, year, isData, sampleType):

    # FOs and padded FOs
    FOs = events.l_fo_conept_sorted
    padded_FOs = ak.pad_none(FOs,3)

    # Filters and cleanups
    filter_flags = events.Flag
    filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.HBHENoiseFilter & filter_flags.HBHENoiseIsoFilter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & ((year == "2016") | filter_flags.ecalBadCalibFilter) & (isData | filter_flags.eeBadScFilter)
    cleanup=events.minMllAFAS > 12

    # IDs
    eleID1=(abs(padded_FOs[:,0].pdgId)!=11) | ((padded_FOs[:,0].convVeto != 0) & (padded_FOs[:,0].lostHits==0))
    eleID2=(abs(padded_FOs[:,1].pdgId)!=11) | ((padded_FOs[:,1].convVeto != 0) & (padded_FOs[:,1].lostHits==0))
    eleID3=(abs(padded_FOs[:,2].pdgId)!=11) | ((padded_FOs[:,2].convVeto != 0) & (padded_FOs[:,2].lostHits==0))

    # Pt requirements for 3rd lepton (different for e and m)
    pt3lmask = ak.any(ak.where(abs(FOs[:,2:3].pdgId)==11,FOs[:,2:3].conept>15.0,FOs[:,2:3].conept>10.0),axis=1)

    # 3l requirements:
    trilep = (ak.num(FOs)) >=3
    pt251510 = (ak.any(FOs[:,0:1].conept > 25.0, axis=1) & ak.any(FOs[:,1:2].conept > 15.0, axis=1) & pt3lmask)
    exclusive = ak.num( FOs[FOs.isTightLep],axis=-1)<4
    mask = (filters & cleanup & trilep & pt251510 & exclusive & eleID1 & eleID2 & eleID3 )

    # MC matching requirement (already passed for data)
    if sampleType == 'prompt':
        lep1_match=((padded_FOs[:,0].genPartFlav==1) | (padded_FOs[:,0].genPartFlav == 15))    
        lep2_match=((padded_FOs[:,1].genPartFlav==1) | (padded_FOs[:,1].genPartFlav == 15))
        lep3_match=((padded_FOs[:,2].genPartFlav==1) | (padded_FOs[:,2].genPartFlav == 15))
        mask = mask & lep1_match & lep2_match & lep3_match
    elif sampleType =='conversions':
        lep1_match=(padded_FOs[:,0].genPartFlav==22)
        lep2_match=(padded_FOs[:,1].genPartFlav==22)
        lep3_match=(padded_FOs[:,2].genPartFlav==22)
        mask = mask & ( lep1_match | lep2_match | lep3_match ) 
    elif sampleType == 'nonprompt':
        lep1_match=((padded_FOs[:,0].genPartFlav!=1) & (padded_FOs[:,0].genPartFlav != 15) & (padded_FOs[:,0].genPartFlav != 22))
        lep2_match=((padded_FOs[:,1].genPartFlav!=1) & (padded_FOs[:,1].genPartFlav != 15) & (padded_FOs[:,1].genPartFlav != 22))
        lep3_match=((padded_FOs[:,2].genPartFlav!=1) & (padded_FOs[:,2].genPartFlav != 15) & (padded_FOs[:,2].genPartFlav != 22))
        mask = mask & ( lep1_match | lep2_match | lep3_match ) 
    elif sampleType == "data":
        pass
    else:
        raise Exception(f"Error: Unknown sampleType {sampleType}.")

    events['is3l'] = ak.fill_none(mask,False)

    # SFs
    events['sf_3l'] = padded_FOs[:,0].sf_nom*padded_FOs[:,1].sf_nom*padded_FOs[:,2].sf_nom
    events['sf_3l_hi'] = padded_FOs[:,0].sf_hi*padded_FOs[:,1].sf_hi*padded_FOs[:,2].sf_hi
    events['sf_3l_lo'] = padded_FOs[:,0].sf_lo*padded_FOs[:,1].sf_lo*padded_FOs[:,2].sf_lo

    # SR:
    events['is3l_SR'] = (padded_FOs[:,0].isTightLep)  & (padded_FOs[:,1].isTightLep) & (padded_FOs[:,2].isTightLep)
    events['is3l_SR'] = ak.fill_none(events['is3l_SR'],False)


# 4l selection
def add4lMaskAndSFs(events, year, isData):

    # FOs and padded FOs
    FOs = events.l_fo_conept_sorted
    padded_FOs = ak.pad_none(FOs,4)

    # Filters and cleanups
    filter_flags = events.Flag
    filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.HBHENoiseFilter & filter_flags.HBHENoiseIsoFilter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & ((year == "2016") | filter_flags.ecalBadCalibFilter) & (isData | filter_flags.eeBadScFilter)
    cleanup = events.minMllAFAS > 12

    # IDs
    eleID1 = ((abs(padded_FOs[:,0].pdgId)!=11) | ((padded_FOs[:,0].convVeto != 0) & (padded_FOs[:,0].lostHits==0)))
    eleID2 = ((abs(padded_FOs[:,1].pdgId)!=11) | ((padded_FOs[:,1].convVeto != 0) & (padded_FOs[:,1].lostHits==0)))
    eleID3 = ((abs(padded_FOs[:,2].pdgId)!=11) | ((padded_FOs[:,2].convVeto != 0) & (padded_FOs[:,2].lostHits==0)))
    eleID4 = ((abs(padded_FOs[:,3].pdgId)!=11) | ((padded_FOs[:,3].convVeto != 0) & (padded_FOs[:,3].lostHits==0)))

    # Pt requirements for 3rd and 4th leptons (different for e and m)
    pt3lmask = ak.any(ak.where(abs(FOs[:,2:3].pdgId)==11,FOs[:,2:3].conept>15.0,FOs[:,2:3].conept>10.0),axis=1)
    pt4lmask = ak.any(ak.where(abs(FOs[:,3:4].pdgId)==11,FOs[:,3:4].conept>15.0,FOs[:,3:4].conept>10.0),axis=1)

    # 4l requirements:
    fourlep  = (ak.num(FOs)) >= 4
    pt25151510 = (ak.any(FOs[:,0:1].conept > 25.0, axis=1) & ak.any(FOs[:,1:2].conept > 15.0, axis=1) & pt3lmask & pt4lmask)
    tightleps = ((padded_FOs[:,0].isTightLep) & (padded_FOs[:,1].isTightLep) & (padded_FOs[:,2].isTightLep) & (padded_FOs[:,3].isTightLep))
    mask = (filters & cleanup & fourlep & pt25151510 & tightleps & eleID1 & eleID2 & eleID3 & eleID4)
    events['is4l'] = ak.fill_none(mask,False)

    # SFs:
    events['sf_4l'] = padded_FOs[:,0].sf_nom*padded_FOs[:,1].sf_nom*padded_FOs[:,2].sf_nom*padded_FOs[:,3].sf_nom
    events['sf_4l_hi'] = padded_FOs[:,0].sf_hi*padded_FOs[:,1].sf_hi*padded_FOs[:,2].sf_hi*padded_FOs[:,3].sf_hi
    events['sf_4l_lo'] = padded_FOs[:,0].sf_lo*padded_FOs[:,1].sf_lo*padded_FOs[:,2].sf_lo*padded_FOs[:,3].sf_lo

    # SR: Don't really need this for 4l, but define it so we can treat 4l category similar to 2lss and 3l
    events['is4l_SR'] = tightleps
    events['is4l_SR'] = ak.fill_none(events['is4l_SR'],False)


def addLepCatMasks(events):

    # FOs and padded FOs
    fo = events.l_fo_conept_sorted
    padded_fo = ak.pad_none(fo,4)
    padded_fo_id = padded_fo.pdgId

    # Find the numbers of e and m in the event
    is_e_mask = (abs(padded_fo_id)==11)
    is_m_mask = (abs(padded_fo_id)==13)
    n_e_2l = ak.sum(is_e_mask[:,0:2],axis=-1) # Make sure we only look at first two leps
    n_m_2l = ak.sum(is_m_mask[:,0:2],axis=-1) # Make sure we only look at first two leps
    n_e_3l = ak.sum(is_e_mask[:,0:3],axis=-1) # Make sure we only look at first three leps
    n_m_3l = ak.sum(is_m_mask[:,0:3],axis=-1) # Make sure we only look at first three leps
    n_e_4l = ak.sum(is_e_mask,axis=-1)        # Look at all the leps
    n_m_4l = ak.sum(is_m_mask,axis=-1)        # Look at all the leps

    # 2l masks
    events['is_ee'] = ((n_e_2l==2) & (n_m_2l==0)) 
    events['is_em'] = ((n_e_2l==1) & (n_m_2l==1)) 
    events['is_mm'] = ((n_e_2l==0) & (n_m_2l==2)) 

    # 3l masks
    events['is_eee'] = ((n_e_3l==3) & (n_m_3l==0)) 
    events['is_eem'] = ((n_e_3l==2) & (n_m_3l==1)) 
    events['is_emm'] = ((n_e_3l==1) & (n_m_3l==2)) 
    events['is_mmm'] = ((n_e_3l==0) & (n_m_3l==3)) 

    # 4l masks
    events['is_eeee'] = ((n_e_4l==4) & (n_m_4l==0))
    events['is_eeem'] = ((n_e_4l==3) & (n_m_4l==1))
    events['is_eemm'] = ((n_e_4l==2) & (n_m_4l==2))
    events['is_emmm'] = ((n_e_4l==1) & (n_m_4l==3))
    events['is_mmmm'] = ((n_e_4l==0) & (n_m_4l==4))
    events['is_gr4l'] = ((n_e_4l+n_m_4l)>4)


# Returns a mask for events with a same flavor opposite sign pair close to the Z
# Mask will be True if any combination of 2 leptons from within the given collection satisfies the requirement
def get_Z_peak_mask(lep_collection,pt_window):
    ll_pairs = ak.combinations(lep_collection, 2, fields=["l0","l1"])
    zpeak_mask = (abs((ll_pairs.l0+ll_pairs.l1).mass - 91.2)<pt_window)
    sfos_mask = (ll_pairs.l0.pdgId == -ll_pairs.l1.pdgId)
    sfosz_mask = ak.flatten(ak.any((zpeak_mask & sfos_mask),axis=1,keepdims=True)) # Use flatten here because it is too nested (i.e. it looks like this [[T],[F],[T],...], and want this [T,F,T,...]))
    return sfosz_mask

def GetJetVariables(jets):
  jj_pairs = ak.combinations(jets, 2, fields=["j0","j1"])
  bjets = jets[(jets.isBtag)]
  
  #ujets = jets[(jets.isBtag==0)]
  #uu_pairs = ak.combinations(ujets, 2, fields=["j0", "j1"])
  mjj = (jj_pairs.j0+jj_pairs.j1).mass
  ptjj = (jj_pairs.j0+jj_pairs.j1).pt
  dR = jj_pairs.j0.delta_r(jj_pairs.j1)

  # Median DR between jets
  dRsort = dR[ak.argsort(dR, axis=-1,ascending=False)]
  medianIndices = np.array(np.ceil( (ak.num(dRsort)-1)/2 ), dtype=int)
  allIndices = ak.local_index(dRsort)
  indexmask = (allIndices == medianIndices)
  dRmedian = dRsort[indexmask]

  # vars for jets with min DR
  args = ak.singletons(ak.argmin(dR, axis=1))
  dR = dR[args]
  mjjmindr = (mjj[args])
  ptjjmindr = ptjj[args]
  j0 = jets[ak.argmax(jets.pt,axis=-1,keepdims=True)]

  return j0, dR, dRmedian, mjjmindr, ptjjmindr

def GetMaxBJet(jets):
    b1=jets[ak.argmax(jets.btagDeepB,axis=-1,keepdims=True)]
    return b1


def GetJetLepVar(jets, leps):
  btags = jets[jets.isBtag]
  jvec = ak.with_name(jets , 'PtEtaPhiMLorentzVector').pvec
  lvec = ak.with_name(leps , 'PtEtaPhiMLorentzVector').pvec
  bvec = ak.with_name(btags, 'PtEtaPhiMLorentzVector').pvec
  jetsum = ak.sum(jvec, axis=-1)
  lepsum = ak.sum(lvec, axis=-1)
  vall = ak.concatenate([jvec, lvec], axis=1)
  vlbb = ak.concatenate([bvec, lvec], axis=1)
  sumall = ak.sum(vall, axis=-1)
  sumlbb = ak.sum(vlbb, axis=-1)
  getpt = lambda v : np.sqrt(v.x*v.x + v.y*v.y)
  ptSumVecAll = getpt(sumall)
  ptSumVeclb  = getpt(sumlbb)
  st = ak.sum(jets.pt, axis=-1) + ak.sum(leps.pt, axis=-1)
  btags['charge'] = np.zeros_like(btags.pt, dtype=int) # charge = 0 just to distinguish from leptons
  lb = ak.concatenate([leps, btags], axis=1)
  lbpairs = ak.combinations(lb, 2, fields=["p1","p2"])
  blmask = np.abs(lbpairs.p1.charge+lbpairs.p2.charge) ==1 # Good pairs formed by btags + leps
  lbpairs = lbpairs[blmask]
  dphi = lambda a,b : (a - b + np.pi) % (2 * np.pi) - np.pi # https://github.com/CoffeaTeam/coffea/blob/7dd4f863837a6319579f078c9e445c61d9106943/coffea/nanoevents/methods/vector.py#L66
  dr = lambda p : np.hypot( p.p1.eta - p.p2.eta, dphi(p.p1.phi, p.p2.phi))
  dRlb = dr(lbpairs)
  dRlb = dRlb[ak.argmin(dRlb,axis=-1,keepdims=True)]
  return ptSumVecAll, ptSumVeclb, dRlb, st

def GetMT(lepton, met):
    ''' Transverse mass with met and lepton '''
    met_pt = (ak.ones_like(lepton.pt))*(met.pt)
    met_px = (ak.ones_like(lepton.pt))*(met.px)
    met_py = (ak.ones_like(lepton.pt))*(met.py)

    return np.sqrt( np.square(lepton.pt + met_pt) - np.square(lepton.px + met_px) - np.square(lepton.py + met_py) )



def GetMlb0b(lepton, jets):
   ''' Invariant mass of lepton and b -- WARNING: run this after applying cuts '''
   l0 = ak.flatten(lepton)
   j0 = jets[ak.argmax(jets.pt,axis=-1,keepdims=True)]
   j0 = ak.flatten(j0)
   mlb = (j0+l0).mass
   return mlb

def GetMlb1b(lepton, jets):
   ''' Invariant mass of lepton and b -- WARNING: run this after applying cuts '''
   # 1b --> Only one posibility
   jet_b = jets[(jets.isBtag)]
   jet_b = ak.flatten(jet_b)
   lep = ak.flatten(lepton)
   if len(jet_b) != len(lep):
     print("WARNING: Number of bjets and leptons are not the same")
     print('Number of bjets: ', len(jet_b), ' .. Number of leptons: ', len(lep), '\n')
   mlb = (jet_b+lep).mass
   return mlb 

def GetMlb2b(lepton, jets):
   ''' Invariant mass of lepton and b -- WARNING: run this after applying cuts '''
   jets_b = jets[jets.isBtag]
   l, b = ak.unzip(ak.cartesian([lepton, jets_b], axis=1))
   mlb_2b = (l+b).mass
   argmin = ak.unflatten(ak.argmin(mlb_2b, axis=1), np.ones_like(ak.argmin(mlb_2b, axis=1)))
   mlb_2b = mlb_2b[argmin]
   mlb = ak.flatten(mlb_2b)
   return mlb

def GetMlb(lepton, jets, nb):
    ''' WARNING: run this after applying cuts '''
    if nb == 0:
        mlb = GetMlb0b(lepton, jets)
    elif nb == 1:
        mlb = GetMlb1b(lepton, jets)
    elif nb == 2:
        mlb = GetMlb2b(lepton, jets)
    return mlb

#def GetMlb(lepton, jets, btagvar='isBtag', nbs=0):
''' Invariant mass of lepton and b '''
'''
   # Masks
   nbtags = ak.num(jets[jets[btagvar]])
   mask_0b = (nbtags == 0)
   mask_1b = (nbtags == 1)
   mask_2b = (nbtags >= 1)

   # 0b --> take leading jet
   if nbs == 0:
     jets_0b = jets[mask_0b]
     leps_0b = lepton[mask_0b]
     j0 = jets_0b[ak.argmax(jets_0b.pt,axis=-1,keepdims=True)]
     l0 = ak.flatten(leps_0b)
     j0 = ak.flatten(j0)
     mlb_0b = (j0+leps_0b).mass
     mlbi = ak.flatten(mlb_0b)

   # 1b --> Only one posibility
   elif nbs ==8:
     jets_1b = jets[mask_1b]
     leps_1b = lepton[mask_1b]
     jets_1b = jets_1b[(jets_1b.isBtag)]
     jet_b = ak.flatten(jets_1b)
     lepts_1b = ak.flatten(leps_1b)
     mlb_1b = (jet_b+leps_1b).mass
     mlbi = ak.flatten(mlb_1b)

   # 2b --> Calculate all the combinations and get the one with minimum m(l,b)
   else:
     jets_g2b = jets[mask_2b]
     leps_g2b = lepton[mask_2b]
     jets_b = jets_g2b[jets_g2b.isBtag]
     l, b = ak.unzip(ak.cartesian([leps_g2b, jets_b], axis=1))
     mlb_2b = (l+b).mass
     argmin = ak.unflatten(ak.argmin(mlb_2b, axis=1), np.ones_like(ak.argmin(mlb_2b, axis=1)))
     mlb_2b = mlb_2b[argmin]
     mlbi = ak.flatten(mlb_2b)

   return mlbi
'''
