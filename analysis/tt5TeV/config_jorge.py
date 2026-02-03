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
from cafea.modules.GetValuesFromJsons import get_lumi
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


# Convert string to list
if   isinstance(ch, str) and ',' in ch: ch = ch.replace(' ', '').split(',')
elif isinstance(ch, str): ch = [ch]
#lumi =  get_lumi(muestras(1))*1000. # pb


year ='2022'
baseweb = 'test/'

#lumi=8*1000 #(fb)
#lumi =  get_lumi(year)*1000. # pb
lumi=8000
#lumi=1
#print(lumi)
#lumi=0
processDic = {
  'tchan': 'TbarBQ_t-channel, TBbarQ_t-channel',
  'tt': 'TTto2L2Nu_TuneCP5,TTtoLNu2Q_TuneCP5',
  #'tW': 'tbarW, tW',
  'tW' : 'TbarWplusto2L2Nu_TuneCP5, TWminusto2L2Nu_TuneCP5, TWminustoLNu2Q_TuneCP5, TbarWplustoLNu2Q_TuneCP5',
  
  #'tW': 'TbarWplusto2L2Nu_TuneCP5',
  #'Non-prompt':'TTToSemiLeptonic, WJetsToLNu_MLM',
  #'DY': 'DYJetsToLL_M50, DYJetsToLL_M10to50', 
  'DY' : 'DYto2L-2Jets_MLL-10to50, DYto2L-2Jets_MLL-50',
  'WJets' : 'WtoLNu-2Jets_TuneCP5',
  'QCD': 'QCD',
  #'QCD_sim' : 'QCD_PT-1000_MuEnrichedPt5,QCD_PT-120to170_MuEnrichedPt5,QCD_PT-15to20_MuEnrichedPt5,QCD_PT-170to300_MuEnrichedPt5,QCD_PT-20to30_MuEnrichedPt5,QCD_PT-300to470_MuEnrichedPt5,QCD_PT-30to50_MuEnrichedPt5,QCD_PT-470to600_MuEnrichedPt5,QCD_PT-50to80_MuEnrichedPt5,QCD_PT-600to800_MuEnrichedPt5,QCD_PT-800to1000_MuEnrichedPt5,QCD_PT-80to120_MuEnrichedPt5, QCD_PT-10to30_EMEnriched,QCD_PT-120to170_EMEnriched,QCD_PT-170to300_EMEnriched,QCD_PT-300toInf_EMEnriched,QCD_PT-30to50_EMEnriched,QCD_PT-50to80_EMEnriched,QCD_PT-80to120_EMEnriched,QCD_PT-15to20_bcToE, QCD_PT-170to250_bcToE,QCD_PT-20to30_bcToE,QCD_PT-250toInf_bcToE,QCD_PT-30to80_bcToE,QCD_PT-80to170_bcToE',
  #'Diboson' : 'WW, WZ, ZZ',#'WWTo2L2Nu, WZTo3LNu',#'WW, WZ, ZZTo2L2Nu',
  #'VV' : 'ZZto2L2Nu_TuneCP5,ZZto2L2Q_TuneCP5, ZZto4L_TuneCP5,ZZZ_TuneCP5,WZZ_TuneCP5, WZto3LNu_TuneCP5,WZto2L2Q_TuneCP5,WWZ_4F, WWW_4F, WWtoLNu2Q_TuneCP5,WWto2L2Nu_TuneCP5 ',#'WWTo2L2Nu, WZTo3LNu',#'WW, WZ, ZZTo2L2Nu',
  'VV' : 'ZZto2L2Nu_TuneCP5,ZZto2L2Q_TuneCP5, ZZto4L_TuneCP5,ZZZ_TuneCP5,WZZ_TuneCP5, WZto3LNu_TuneCP5,WZto2L2Q_TuneCP5,WWZ_4F, WWW_4F, WWtoLNu2Q_TuneCP5,WWto2L2Nu_TuneCP5,WZtoLNu2Q_TuneCP5',#'WW, WZ, ZZTo2L2Nu', #all samples
  #'VV' : 'ZZto2L2Nu_TuneCP5',
  'ttX' : 'TTG-1Jets_PTG-100to200,TTG-1Jets_PTG-10to100,TTG-1Jets_PTG-200,TTZ-ZtoQQ-1Jets_TuneCP5',
  #'VV': 'WZto3LNu_TuneCP5,WZto2L2Q_TuneCP5,ZZZ_TuneCP5, WWto2L2Nu_TuneCP5,WZto2L2Q_TuneCP5 ',
  #'data' : 'MuonEG,EGamma,DoubleMuon,SingleMuon,Muon'
  #'data' : 'MuonEG, DoubleMuon, SingleMuon, SingleElectron, DoubleEG, EGamma',
  'data' : 'MuonEG, DoubleMuon, SingleMuon, EGamma, Muon',
  #'data': 'MuonEG',
}

processDic_noQCD = processDic.copy()
processDic_noQCD.pop('QCD')

bkglist = ['tchan','tt', 'tW', 'DY', 'WJets','VV','ttX','QCD']#, 'QCD_sim']#,'ttX' 'QCD'] #
bkglist_noQCD = ['tchan','tt', 'WJets','tW', 'DY']#, 'VV', 'ttX']#, 'QCD_sim']#,'ttX']#'tchan', 
#bkglist_noQCD = ['QCD_sim']
#bkglist    = ['tt', 'tW', 'Non-prompt','DY', 'VV','ttX']
#bkglist = list(processDic.keys())
bkgnormunc = [0.05, 0.10, 0.3, 0.2, 0.3, 0.3]#,0.3, 0.5] 
diclegendlabels = {'None':'Data', 'tt':'$\\mathrm{t\\bar{t}}$', 'DY':'Drell-Yan','tchan':r'$t$ channel', 'WJetsL':'W+jets (l)','WJetsH':'W+jets (h)', 'QCD':'QCD'}
'''
colordic ={
  'tchan' : '#FFA90E',
  'tt' : '#cc0000',
  'tW' : '#ffc207',
  'semilep': '#6c3b2a',
  #'WJets': '#47ce33', #verde
  'WJets': '#3F90DA',
  'DY': '#3b78cb',
  'VV' : '#cc6235',
  #'ttX' : '#cc6235',
  }
  
'''
colordic ={
  'tt' : '#BD1F01',
  'tW' : '#A96B59', 
  'tchan' :'#FFA90E',
  'WJets': '#92DADD',
  'VV':'#3F90DA',  
  'DY':  '#94A4A2',
  'ttX': '#E76300',
  'QCD' : '#717581',
  'QCD_sim': '#717581'
}



colors = [colordic[k] for k in colordic.keys()]
def GetChLab(channel):
  if '_fake' in channel: 
    channel = 'non-iso ' + channel[0]
  channel = channel.replace('m', '$\mu$')
  return channel

def GetLevLab(lev):
  if   lev == 'dilep'  : return ''
  elif lev == 'g2jets': return '$\geq$2 jets'
  elif lev == 'g2jetsg1b': return ', $\geq$2 jets, $\geq$1 b jet'
  return lev

def GetModSystHisto(path, fname, systname, var=None, prname='tt', samplab='sample', prlab='process', systlab='syst', systnormlab='norm'):
  h  = GetHisto(path+ fname +   '.pkl.gz', var, group=None)
  axes = [x.name for x in h.sparse_axes()]
  if not samplab in axes: return
  sampName = h.identifiers(samplab)[0]
  h = GroupKeepOrder(h, [[systlab, systlab, {systname:systnormlab}], [samplab, prlab, {prname:sampName}] ])
  return h

def GetModSystHistos(path, fname, systname, var=None):
  if fname=='TTTo2L2Nu_hdamp':
    up = GetModSystHisto(path, fname+'UP'+'_'+year,   systname+'Up', var)
    do = GetModSystHisto(path, fname+'DOWN_'+year, systname+'Down', var)
  elif fname=='TTTo2L2Nu_mtop':
    up = GetModSystHisto(path, fname+'175p5'+'_'+year,   systname+'Up', var)
    do = GetModSystHisto(path, fname+'169p5'+'_'+year, systname+'Down', var)
  elif fname=='TTTo2L2Nu_TuneCP5':
    up = GetModSystHisto(path, fname+'up'+'_'+year,   systname+'Up', var)
    do = GetModSystHisto(path, fname+'down'+'_'+year, systname+'Down', var)
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
  #elif var == "njets" and 'level' != 'incl':
    #b0 = 4; bN = 10
  elif var in ['st']:
    b0 = 120; bN = 600;
  elif var in ['sumallpt']:
    b0 = 0; bN = 200
    xtit = '$\sum_\mathrm{j,\ell}\,\mathrm{p}_{T}$ (GeV)'
  elif var in ['u0pt', 'ptuu', 'ptjj']:
    b0 = 2;
  #elif var in ['metnocut']:
    #b0 = 4;
  elif var in ['MVAscore']:
    b0 = 0.2; bN = 0.8
    binRebin = 2
  elif var in ['ht']:
    b0 = 100; bN = 450
    binRebin = 2;
  #elif var in ['j0pt']:
    #b0 = 40; bN = 200
  elif var in ['mjj', 'muu']:
    b0 = 25; bN = 150
  #elif var in ['mlb']:
    #b0 = 25; bN = 200
  elif var in ['dRlb']:
    b0 = 0.5; bN = 2.9
  #elif var in ['njetsnbtags12']:
    #b0 =9 #; bN = 16
    #binRebin = 10
  if b0 is not None:
    p.SetRebin(var, b0, bN, includeLower=True, includeUpper=True, binRebin=binRebin)
  return xtit
