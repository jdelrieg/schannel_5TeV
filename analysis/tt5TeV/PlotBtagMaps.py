from cafea.modules.paths import cafea_path
from coffea import hist, lookup_tools
from cafea.plotter.plotter import GetHisto, GetSFfromCountsHisto, DrawEff, DrawEff2D, GetH2DfromXY, loadHistos, GetEff
import matplotlib.pyplot as plt
import numpy as np
from config import *


def DrawHistoBtag(path, WP='medium', flav=5, year=None, outpath=''):
  var='jetptetaflav'
  path = path%year if year is not None else path
  hu = GetHisto(path, var, {'WP':WP,   'flav':flav})#  , 'flav':flav})
  hd = GetHisto(path, var, {'WP':'all', 'flav':flav})
  hnum = GetHisto(path, var, {'WP':WP,})
  hden = GetHisto(path, var, {'WP':'all'})
  X, ratio = GetEff(hu, hd)
  ratio = ratio[0]
  cy = hu.axis('pt').edges()
  cx = hu.axis('abseta').edges()
  #if year in ['2017', '2018']: ratio = ratio[1:,:-1]
  #else: ratio = ratio[:,:-1]
  #ratio = ratio[1:,:-1]
  ratio = ratio[0:,:-1] #Removing overflow in eta (we dont cover above than 2.5)
  ratio=np.delete(ratio,-1,axis=0) #Removing overflowin pt
  h2d = GetH2DfromXY([cx, cy], ratio, ytit='Y', xtit=['X1', 'X2'], hname=['pt','eta'])
  getnum = lookup_tools.dense_lookup.dense_lookup(hnum.values(overflow='over')[()], [hnum.axis('pt').edges(), hnum.axis('abseta').edges(), hnum.axis('flav').edges()])
  getden = lookup_tools.dense_lookup.dense_lookup(hden.values(overflow='over')[()], [hden.axis('pt').edges(),  hnum.axis('abseta').edges(),hden.axis('flav').edges()])
  pt  = np.array([25., 40., 50.]);
  eta = np.array([0.1, 0.1, 0.1]);
  fla = np.array([1, 4, 5]);
  sflav = 'b' if flav==5 else ('c' if flav==4 else 'l')
  DrawEff2D(h2d,'pt', error=None, error2=None, ytit='$p_{T}$ (GeV)', xtit='|$\eta$|', tit='', outname=outpath+('btagSF_%s_%s_%s.png'%(year, WP, sflav) if year is not None else 'btagSF_%s_%s.png'%(WP, sflav)))
  DrawEff2D(h2d,'pt', error=None, error2=None, ytit='$p_{T}$ (GeV)', xtit='|$\eta$|', tit='', outname=outpath+('btagSF_%s_%s_%s.pdf'%(year, WP, sflav) if year is not None else 'btagSF_%s_%s.pdf'%(WP, sflav)))


outpath = baseweb+datatoday+'/B-tagging/'
outpath='btag_2p5_onebin/'
if not os.path.exists(outpath):
  os.makedirs(outpath)

# cafea/data/btagSF/UL/btagMCeff_5TeV.pkl.gz
if not os.path.isfile(path):
  print('File not found. Specify SF files with the -p options.')
  exit()

for wp in ['loose', 'medium', 'tight']:
  for f in [1, 4, 5]:
    DrawHistoBtag(path, wp, f, outpath=outpath)

    

    
