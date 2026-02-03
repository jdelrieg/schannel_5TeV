from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import gzip
import pickle
import json
import os
import uproot
import matplotlib.pyplot as plt
import numpy as np
from coffea import hist, processor 
from coffea.hist import plot
import os, sys

from cafea.plotter.plotter import plotter

import argparse
parser = argparse.ArgumentParser(description='You can customize your run')
parser.add_argument('--path',     '-p', default = 'histos/plots5TeV.pkl.gz', help = 'Path to pkl file')
parser.add_argument('--variable', '-v', default = None                 , help = 'Variable')
parser.add_argument('--channel',  '-c', default = 'em'                     , help = 'Channels')
parser.add_argument('--level',    '-l', default = 'incl'                   , help = 'Variable')
parser.add_argument('--output',   '-o', default = None                     , help = 'Name of the output png file')
parser.add_argument('--data',     '-d', action= 'store_true'             , help = 'Do data?')
args = parser.parse_args()

path  = args.path
var = args.variable
ch = args.channel
level = args.level
output = args.output
doData = args.data
syst = 'norm'

# Convert string to list
if   isinstance(ch, str) and ',' in ch: ch = ch.replace(' ', '').split(',')
elif isinstance(ch, str): ch = [ch]

lumi = 302; # pb
year = 2017

processDic = {
  'tt': 'tt',
  'tW': 'tbarW, tW',
  'WJets': 'W0JetsToLNu, W1JetsToLNu, W2JetsToLNu, W3JetsToLNu',
  'DY': 'DYJetsToLLMLL50, DYJetsToLLM10to50',
  'data' : 'SingleMuon, HighEGJet',
}
bkglist = ['tt', 'tW', 'WJets', 'DY']

colordic ={
  'tt' : '#ec0000',
  'tW' : '#ff9c00',
  'WJets': '#00a714',
  'DY': '#8f8f8f',
  'QCD' : '#0096ff',
}

colors = [colordic[k] for k in bkglist]

def GetQCDbkg(var, categories):
  ch = categories['channel']
  l  = categories['level']
  plt = plotter(path, prDic=processDic, bkgList=bkglist, lumi=lumi)
  cat_fake = categories.copy()
  cat_fake['channel'] = [(c+'_fake' if not 'fake' in c else c) for c in cat_fake['channel']]
  # Non isolated data - prompt MC bkg
  h_data_fake = plt.GetHistogram(var, ['data'], cat_fake).group('process', hist.Cat("process", "process"), {'QCD':'data'})
  h_mc_fake   = plt.GetHistogram(var, bkglist, cat_fake).group('process', hist.Cat("process", "process"), {'QCD':bkglist})
  h_mc_fake.scale(-1*lumi)
  htot = (h_data_fake + h_mc_fake)

  # Factor using MET < 20 GeV
  data_metl20    = plt.GetYields('counts_metl20', categories, pr='data')
  mc_metl20      = plt.GetYields('counts_metl20', categories, pr=bkglist)
  data_metl20_fk = plt.GetYields('counts_metl20', cat_fake,   pr='data')
  mc_metl20_fk   = plt.GetYields('counts_metl20', cat_fake,   pr=bkglist)
  fact = (data_metl20 - mc_metl20)/(data_metl20_fk - mc_metl20_fk)
  htot.scale(fact)

  # Remove negative values
  #w = [(-1*x if x < 0 else 0.) for x in htot.values()[list(htot.values().keys())[0]] ]
  #centers = htot.axis(var).centers()
  #htot.fill(met=np.array(centers), process='QCD', weight=np.array(w))

  return htot

def Draw(var, categories, output=None, label='', outpath='temp/', doQCD=False):
  plt = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi)
  y = plt.GetYields('counts', pr='data')
  h = plt.hists['met']
  if doQCD:
    hQCD = GetQCDbkg(var, categories)
    plt.AddExtraBkgHist(hQCD)
  plt.SetOutpath(outpath)
  plt.plotData = doData
  plt.SetCategories(categories)
  plt.SetRegion(label)
  plt.SetOutput(output)
  plt.Stack(var, xtit='', ytit='')
  #plt.PrintYields('counts')

def Print2lplots():
  for c in ['em', 'ee', 'mm']:
    for l in ['incl', 'g2jets']:
      outp = outpath+'/2l/'+c+'/'+l+'/'
      cat = {'channel':c, 'level':l, 'syst':'norm'}
      for var in ['counts', 'njets', 'nbtags', 'ht', 'met', 'j0pt', 'j0eta', 'l0pt', 'l0eta', 'invmass', 'invmass2']:
        if l=='incl' and var in ['j0pt', 'j0eta']: continue
        outname = "%s_%s_%s"%(var, c, l)
        Draw(var, cat, outname, outpath=outp)

def Print1lplots():
  outp = outpath+'/1l/'
  for c in ['e', 'm', 'e_fake', 'm_fake']:
    for l in ['incl', 'g4jets', '0b', '1b', '2b']:
      cat = {'channel':c, 'level':l, 'syst':'norm'}
      outp = outpath+'/1l/'+c+'/'+l+'/'
      for var in ['mlb']:#['mjj', 'mt', 'ptjj', 'minDRjj']: #['counts', 'njets', 'nbtags', 'ht', 'met', 'j0pt', 'j0eta', 'ept', 'eeta', 'mpt', 'meta']:
        if l=='incl' and var in ['j0pt', 'j0eta']: continue
        outname = "%s_%s_%s"%(var, c, l)
        Draw(var, cat, outname, outpath=outp, doQCD=True)


if not var is None:
  categories = { 'channel': ch, 'level' : level, 'syst':'norm'}
  Draw(var, categories, output, doQCD=True)


else:
  outpath = '/nfs/fanae/user/jriego/'
  doData = True
  #Print2lplots()
  Print1lplots()


