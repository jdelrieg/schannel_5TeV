#!/usr/bin/env python3
from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import gzip
import pickle
import json
import uproot3
import matplotlib.pyplot as plt
import numpy as np
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea import hist, processor 
from coffea.hist import plot
import os, sys
from cafea.plotter.OutText import OutText

from cafea.plotter.plotter import plotter, GetH1DfromXY
from cafea.plotter.plotter import *


sys.path.append('/nfs/fanae/user/jriego/.conda/envs/conda-jriego-aug23/lib/python3.9/site-packages')

import argparse
parser = argparse.ArgumentParser(description='You can customize your run')
parser.add_argument('--path',     '-p', default = 'histos/plots5TeV.pkl.gz', help = 'Path to pkl file')
parser.add_argument('--variable', '-v', default = None                 , help = 'Variable')
parser.add_argument('--channel',  '-c', default = 'em'                     , help = 'Channels')
parser.add_argument('--level',    '-l', default = 'incl'                   , help = 'Variable')
parser.add_argument('--output',   default = None                     , help = 'Name of the output png file')
parser.add_argument('--outpath',  '-o', default = None                     , help = 'Name of the output path')
parser.add_argument('--data',     '-d', action= 'store_true'             , help = 'Do data?')
parser.add_argument('--syst',     '-s', default= None             , help = 'Systematic choice')
parser.add_argument('--nSlots',   '-n', default= 4             , help = 'Number of slots for parallelization')
parser.add_argument('--verbose',   default= 0             , help = 'level of verbosity')
parser.add_argument('--force',  '-f', action= 'store_true'             , help = 'Force to overwrite')
parser.add_argument('--inputFile',  default=''             , help = 'Used for combine scripts')
args = parser.parse_args()

path  = args.path
var = args.variable
ch = args.channel
level = args.level
output = args.output
doData = args.data
outpatho = args.outpath
systch = args.syst
verbose = int(args.verbose)
nSlots = int(args.nSlots)
inputFile = args.inputFile
force = args.force
if outpatho is None: outpatho = 'temp/'
if not outpatho.endswith('/'): outpatho += '/'
#syst = 'norm'

pathQCD = path

from datetime import datetime
now = datetime.now()
datatoday = str(now.strftime('%d')) + str(now.strftime('%B')).lower()[:3] + str(now.strftime('%Y'))[2:]
datatoday = '30aug2023_tch'
#datatoday = 'splitJES_0b_UE_met30_xTrigSF_splitPDFs_lumiUnc_mujetsB_JER_METfilters_tchannel'
baseweb = '/nfs/fanae/user/jriego/www/public/tchannel5TeV/'


# Convert string to list
if   isinstance(ch, str) and ',' in ch: ch = ch.replace(' ', '').split(',')
elif isinstance(ch, str): ch = [ch]
if   isinstance(level, str) and ',' in level: level = level.replace(' ', '').split(',')
lumi = 302; # pb
year = 2017

processDic = {
  'tchan': 'tbarchannel, tchannel',
  'tt': 'ttPS',#, ttPS',
  'tW': 'tbarW, tW',
  
  'WJets':  'W0JetsToLNu, W1JetsToLNu, W2JetsToLNu, W3JetsToLNu',# 'WJetsToLNu',#
  'QCD': 'QCD',
  'DY': 'DYJetsToLLMLL50, DYJetsToLLM10to50',
  'data' : 'SingleMuon, HighEGJet',
}

diclegendlabels = {'None':'Data', 'tchan':'t-channel','tt':'$\\mathrm{t\\bar{t}}$', 'DY':'Drell-Yan', 'WJets':'W+jets', 'QCD':'QCD'}

processDic_noQCD = processDic.copy()
processDic_noQCD.pop('QCD')


bkglist    = ['tchan', 'tt', 'tW',  'WJets', 'DY', 'QCD']#'tchan', 
bkglist_noQCD = ['tchan','tt', 'tW',  'WJets', 'DY']#'tchan', 
bkgnormunc = [0.02,0.05, 0.056,  0.2, 0.2, 0.3]



colordic ={
  'tchan' : '#EEC895',#'#fa00fa',
  'tt' : '#cc0000',
  'tW' : '#ffc207',
  
  'WJets': '#47ce33',
  'DY': '#3b78cb',
  'QCD' : '#aaaaaa',
}

colors = [colordic[k] for k in bkglist]

def GetChLab(channel):
  if isinstance(channel, list) and len(channel) > 1:
    channel = '$\ell$'
  elif isinstance(channel, list):
    channel = channel[0]
  if '_fake' in channel: 
    channel = 'non-iso ' + channel[0]
  channel = channel.replace('m', '$\mu$')
  return channel

def GetLevLab(lev):
  if   lev == 'incl'  : return ''
  elif lev == 'g2jets': return ', $\geq$2 jets'
  elif lev == 'g3jets': return ', $\geq$3 jets'
  elif lev == 'g4jets': return ', $\geq$4 jets'
  elif lev == '0b'    : return ', 0b'
  elif lev == '1b'    : return ', $\geq$4 jets, 1b'
  elif lev == '2b'    : return ', $\geq$4 jets, 2b'
  elif lev == 'g5j1b' : return ', $\geq$5j, 1b'
  elif lev == 'g5j2b' : return ', $\geq$5j, 2b'
  return lev

def GetModSystHisto(path, fname, systname, var=None, prname='tt', samplab='sample', prlab='process', systlab='syst', systnormlab='norm'):
  h  = GetHisto(path+ fname +   '.pkl.gz', var, group=None)
  axes = [x.name for x in h.sparse_axes()]
  if not samplab in axes: return
  sampName = h.identifiers(samplab)[0]
  h = GroupKeepOrder(h, [[systlab, systlab, {systname:systnormlab}], [samplab, prlab, {prname:sampName}] ])
  return h

def GetModSystHistos(path, fname, systname, var=None):
  up = GetModSystHisto(path, fname+'Up',   systname+'Up', var)
  do = GetModSystHisto(path, fname+'Down', systname+'Down', var)
  return up, do


  #elif var in ['ht']:
  #  b0 = 2

def RebinVar(p, var, level=None):
  b0 = None; bN = None; binRebin=None
  xtit = None

  if var in ['minDRjj', 'minDRuu']:
    b0 = 0.4; bN = 2.0
  elif var =='medianDRjj':
    #b0 = 1.0; bN = 3.5
    b0 = 1.4; bN = 3.2
  elif var in ['medianDRuu']:
    b0 = 0.5; bN = 3.7
  elif var == "njets" and 'level' != 'incl':
    b0 = 0; bN = 10
  elif var in ['st']:
    b0 = 120; bN = 600;
  elif var in ['sumallpt']:
    b0 = 0; bN = 200
    xtit = '$\sum_\mathrm{j,\ell}\,\mathrm{p}_{T}$ (GeV)'
  elif var in ['met','u0pt', 'ptuu', 'ptjj']:
    b0 = 2;
  elif var in ['metnocut']:
    b0 = 4;
  elif var in ['MVAscore']:
    b0 = 0.1; bN = 0.9
    #b0 = 0.2; bN = 0.8
    binRebin = 2
  elif var in ['ht']:
    b0 = 100; bN = 450
    binRebin = 2;
  elif var in ['j0pt']:
    b0 = 40; bN = 200
  elif var in ['mjj', 'muu']:
    b0 = 25; bN = 150
  elif var in ['mlb']:
    b0 = 25; bN = 200
  elif var in ['dRlb']:
    b0 = 0.5; bN = 2.9
  if b0 is not None:
    p.SetRebin(var, b0, bN, includeLower=True, includeUpper=True, binRebin=binRebin)
  return xtit
