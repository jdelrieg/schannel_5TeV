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

from cafea.modules.GetValuesFromJsons import get_param
from cafea.analysis.objects import *
from cafea.analysis.corrections import GetBTagSF, GetBtagEff, AttachMuonSF, AttachElectronSF, GetPUSF, GetTriggerSF5TeV, jet_factory, jet_factory_data, met_factory, GetBtagSF5TeV, GetPUSF, AttachMuonPOGSFs, AttachElecPOGSFs, GetTriggerSF, GetTrigSFttbar
from cafea.analysis.selection import *
from cafea.modules.paths import cafea_path
nbins = 1500
lista = list(10 ** np.linspace(np.log10(0.05), np.log10(150), nbins))
ar = np.array(lista)

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples):

        self._samples = samples

        # Create the histograms
        # 'name' : hist.Hist("Ytitle", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat("syst", "syst"), hist.Bin("name", "X axis (GeV)", 20, 0, 100)),
        self._accumulator = processor.dict_accumulator({
        #'invmass'    : hist.Hist("Events", hist.Bin("invmass", "$m_{\mu\mu}$ (GeV) ", 2400, 0, 120)),

        'invmassg'    : hist.Hist("Events", hist.Bin("invmassg", "$m_{\mu\mu}$ (GeV) ", ar)),
        'invmass'    : hist.Hist("Events", hist.Bin("invmass", "$m_{\mu\mu}$ (GeV) ", 500, 0.5, 120)),
        'Z'    : hist.Hist("Events", hist.Bin("Z", "$m_{\mu\mu}$ (GeV) ", 20, 80, 100)),
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
        hout = self.accumulator.identity()

        # Initialize objects
        mu   = events.Muon

        mu['isGood'] = isMuonPOGL(mu, ptCut=0)
        m_sel = mu[mu.isGood] #mu[mu.isLoose & mu.isMVA]
        mmpairs = ak.combinations(m_sel, 2, fields=["m0","m1"])
        mll = (mmpairs.m0+mmpairs.m1).mass # Invmass for leading two leps

        mll_flat = ak.flatten(mll)
        weights = np.ones_like(mll_flat, dtype=bool)

        events['ismm'] = (ak.num(m_sel) == 2)
        hout['invmassg'].fill(invmassg=mll_flat, weight=weights)
        hout['invmass'].fill(invmass=mll_flat, weight=weights)
        hout['Z'].fill(Z=mll_flat, weight=weights)

        return hout
  
    def postprocess(self, accumulator):
        return accumulator

if __name__ == '__main__':
    # Load the .coffea files
    outpath= './coffeaFiles/'
    samples     = load(outpath+'samples.coffea')
    topprocessor = AnalysisProcessor(samples)

