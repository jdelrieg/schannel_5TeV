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
from cafea.analysis.selection import *
from cafea.modules.paths import cafea_path

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples):

        self._samples = samples
        self._accumulator = processor.dict_accumulator({
        'PDF'        : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("PDF",     "Counts", 103, 0, 103)),
        'Scales'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Bin("Scales",  "Counts", 9, 0, 9)),
        'counts'     : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("counts",  "Counts", 1, 0, 10)),
        'njets'      : hist.Hist("Events", hist.Cat("sample", "sample"), hist.Cat("channel", "channel"), hist.Cat("level", "level"), hist.Cat('syst', 'syst'), hist.Cat('sign', 'sign'), hist.Bin("njets",   "Jet multiplicity", 6, 0, 6)),
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
        doPS         = (histAxisName in ['tt', 'ttPS', 'TTTo2L2Nu']) and events.PSWeight is not None and len(events.PSWeight[0])>=4
        doPDFunc = "sumPDFWeights" in self._samples[dataset]

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

        # Get the gen leptons
        lep  = events.GenDressedLepton
        jet  = events.GenJet
        met  = events.GenMET

        # Leptons
        lep = lep[(lep.pt>35) & (np.abs(lep.eta) < 2.4)]
        nleps = ak.num(lep)
        
        # Jets
        jet = jet[(jet.pt>30) & (np.abs(jet.eta) < 2.4)]
        njets = ak.num(jet)

        # Met cut (30 GeV)
        metcut = (met.pt >= 30)

        # Inv mass
        llpairs = ak.combinations(lep, 2, fields=["l0", "l1"])
        mll = (llpairs.l0+llpairs.l1).mass

        mllvalues = np.where(ak.num(mll)==0, [[0]], mll)
        mllvalues = np.where(ak.num(mllvalues)>1, [[0]], mllvalues)
        mllvalues = ak.flatten(mllvalues, axis=1)

        mllcut = (mllvalues>20)
        offZ = (np.abs(mllvalues-90) > 15)

        # Final states
        lpad = ak.pad_none(lep, 2)
        l0 = lpad[:,0]
        l1 = lpad[:,1]
        id0 = l0.pdgId
        id1 = l1.pdgId
        isem = ( (np.abs(id0) == 11) & (np.abs(id1) == 13) ) | ( (np.abs(id0) == 13) & (np.abs(id1) == 11) )
        isee = ( (np.abs(id0) == 11) & (np.abs(id1) == 11) )
        ismm = ( (np.abs(id0) == 13) & (np.abs(id1) == 13) )

        # Weights (ones)
        genw = np.ones_like(events["event"]) #genw = events["genWeight"]
        weights_dict = coffea.analysis_tools.Weights(len(events),storeIndividual=True)
        weights_dict.add("norm", genw)

        if doPS: 
          i_ISRdown = 0; i_FSRdown = 1; i_ISRup = 2; i_FSRup = 3
          ISRUp = ak.flatten(events.PSWeight[ak.Array(ak.local_index(events.PSWeight)==i_ISRup)])
          ISRDo = ak.flatten(events.PSWeight[ak.Array(ak.local_index(events.PSWeight)==i_ISRdown)])
          FSRUp = ak.flatten(events.PSWeight[ak.Array(ak.local_index(events.PSWeight)==i_FSRup)])
          FSRDo = ak.flatten(events.PSWeight[ak.Array(ak.local_index(events.PSWeight)==i_FSRdown)])
          weights_dict.add('ISR', np.ones_like(events["event"]), ISRUp, ISRDo)
          weights_dict.add('FSR', np.ones_like(events["event"]), FSRUp, FSRDo)

        systList = ["norm"]
        if doPS: systList += ['ISRUp', 'ISRDown', 'FSRUp', 'FSRDown']

        counts = np.ones_like(events['event'], dtype=float)

        hout = self.accumulator.identity()
        channels = ['em', 'ee', 'mm'] 
        levels = ['dilep', 'g2jets', 'offZ', 'metcut']

        selections = PackedSelection(dtype='uint64')
        selections.add("em", (isem))
        selections.add("ee", (isee))
        selections.add("mm", (ismm))
        selections.add("dilep",  (nleps >= 2))
        selections.add("g2jets", (njets >= 2))
        selections.add("offZ",   (offZ)&(njets >= 2))
        selections.add("metcut", (njets >= 2)&(offZ)&(metcut))
        selections.add("mll", (mllcut) )

        # Save the total number of events

        weights = weights_dict.weight(None)
        hout['counts'].fill(sample=histAxisName, channel='all', level='all', counts=counts, syst="norm", sign="all", weight=weights)
        for syst in systList:
         for ch in channels:
          for lev in levels:
            cuts = [ch] + [lev] + ['mll', 'dilep']
            cut   = selections.all(*cuts)
            weights = weights_dict.weight(None)
            weights = weights[cut]
            hout['counts'].fill(sample=histAxisName, channel=ch, level=lev, counts=counts[cut],  syst=syst, sign='OS', weight=weights)
            hout['njets'].fill(sample=histAxisName, channel=ch, level=lev, njets=njets[cut], syst=syst, sign='OS', weight=weights)

            # Fill scale and pdf uncertainties
            if doPDFunc:
              scale_w = np.transpose(scaleweights[cut])*(weights)
              pdf_w   = np.transpose(pdfweights  [cut])*(weights)
              hout['Scales'].fill(sample=histAxisName, channel=ch, level=lev, Scales=ak.flatten(scaleweights_bins[cut]), syst="norm", weight=ak.flatten(scale_w))
              hout['PDF']   .fill(sample=histAxisName, channel=ch, level=lev, PDF   =ak.flatten(pdfweights_bins[cut]),   syst="norm", weight=ak.flatten(pdf_w))
        return hout
  
    def postprocess(self, accumulator):
        return accumulator

if __name__ == '__main__':
  processor = AnalysisProcessor()
  print('Execute the proper run.py script!! ')


