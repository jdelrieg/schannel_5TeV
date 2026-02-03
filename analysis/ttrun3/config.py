from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import gzip
import pickle
import json
import uproot3
import matplotlib.pyplot as plt
import numpy as np
from coffea import hist, processor 
from coffea.hist import plot
import os, sys
from cafea.plotter.OutText import OutText

from cafea.plotter.plotter import plotter, GetH1DfromXY
from cafea.plotter.plotter import *

import argparse
parser = argparse.ArgumentParser(description='You can customize your run')
parser.add_argument('--path',     '-p', default = 'histos/plots5TeV.pkl.gz', help = 'Path to pkl file')
parser.add_argument('--variable', '-v', default = None                 , help = 'Variable')
parser.add_argument('--channel',  '-c', default = 'em'                     , help = 'Channels')
parser.add_argument('--level',    '-l', default = 'dilep'                   , help = 'Variable')
parser.add_argument('--output',   default = None                     , help = 'Name of the output png file')
parser.add_argument('--outpath',  '-o', default = None                     , help = 'Name of the output path')
parser.add_argument('--data',     '-d', action= 'store_true'             , help = 'Do data?')
args = parser.parse_args()

path  = args.path
var = args.variable
ch = args.channel
level = args.level
output = args.output
doData = args.data
outpatho = args.outpath
if outpatho is None: outpatho = ''
if not outpatho.endswith('/'): outpatho += '/'
syst = 'norm'

# Convert string to list
if   isinstance(ch, str) and ',' in ch: ch = ch.replace(' ', '').split(',')
elif isinstance(ch, str): ch = [ch]
lumi =  1200#314.6##142.8#93.7 #+142.8 #129 + 93.7; #142.8#129 + 93.7; # pb
year = '2022'


processDic = {
  'tt': 'TTTo2L2Nu',
  'tW': 'tbarW, tW',
  'semilep':'TTToSemiLeptonic',
  'WJets':'WJetsToLNu',
  'DY': 'DYJetsToLL_M50, DYJetsToLL_M10to50', 
  'Diboson' : 'WW, WZ, ZZ',#'WWTo2L2Nu, WZTo3LNu',#'WW, WZ, ZZTo2L2Nu',
  'data' : 'MuonEG,EGamma,DoubleMuon,SingleMuon,Muon'
}

bkglist    = ['tt', 'tW', 'semilep','WJets','DY', 'Diboson']
#bkglist = list(processDic.keys())
bkgnormunc = [0.05, 0.15, 0.2, 0.3, 0.2, 0.3]

colordic ={
  'tt' : '#cc0000',
  'tW' : '#ffc207',
  'semilep': '#6c3b2a',
  'WJets': '#47ce33',
  'DY': '#3b78cb',
  'Diboson' : '#fdffcb',
}


colors = [colordic[k] for k in colordic.keys()]

def GetChLab(channel):
  channel = channel.replace('m', '$\mu$')
  return channel

def GetLevLab(lev):
  if   lev == 'dilep'  : return ''
  elif lev == 'g2jets': return ', $\geq$2 jets'
  elif lev == 'g2jetsg1b': return ', $\geq$2 jets, $\geq$1 b jet'
  return ''

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

def GetJECSystHistos(path, fname, var, categories): # placeholder : up=down
  #up  = GetModSystHisto(path, fname , systname+'Up', var)
  #do  = GetModSystHisto(path, fname , systname+'Down', var).values()
  catecopy =  categories.copy()
  catecopy['sign'] = 'OS'
  nom=GetHisto(path+ 'TTTo2L2Nu.pkl.gz', var, catecopy)
  nom=nom.integrate('syst','norm')
  up=GetHisto(path+ fname+'.pkl.gz', var, catecopy)
  up=up.integrate('syst','norm')
  var=abs(nom.values()[('TTTo2L2Nu',)]-up.values()[('TTTo2L2Nu',)])/nom.values()[('TTTo2L2Nu',)]
  return var


