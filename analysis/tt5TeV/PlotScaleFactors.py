from cafea.modules.paths import cafea_path
from cafea.modules.utils import loadHistos
from coffea import hist, lookup_tools
from cafea.plotter.plotter import GetHisto, GetSFfromCountsHisto, DrawEff, DrawEff2D, GetH2DfromXY
from statsmodels.stats.proportion import proportion_confint
import matplotlib.pyplot as plt
import numpy as np

def GetHistoTriggerSF5TeV(path, ch='m'):
  integrate=['pr']
  dataSample = 'HighEGJet' if ch =='m' else 'SingleMuon'
  hu = GetHisto(path, var, {'channel':ch, 'val':'num'}, group=['sample', 'pr', {'pr':mcSample}], integrate=integrate)
  hd = GetHisto(path, var, {'channel':ch, 'val':'den'}, group=['sample', 'pr', {'pr':mcSample}], integrate=integrate)
  du = GetHisto(path, var, {'channel':ch, 'val':'num'}, group=['sample', 'pr', {'pr':dataSample}], integrate=integrate)
  dd = GetHisto(path, var, {'channel':ch, 'val':'den'}, group=['sample', 'pr', {'pr':dataSample}], integrate=integrate)
  ratio, do, up = GetSFfromCountsHisto(hu, hd, du, dd)
  cx = hu.axis('pt').edges()
  cy = hu.axis('abseta').edges()

def DrawTriggerSFs(path, ch='m', var='pteta', dataSample='HighEGJet', mcSample='tt', xtit='p$_T$ (GeV)', ytit='$|\eta|$', tit='Muon trigger scale factors', outname='temp.png', axesToIntegrate=''):
  ''' Draw trigger SFs for the tt 5 TeV analysis '''
  integrate = ['pr']
  if axesToIntegrate != '':
    axesToIntegrate = axesToIntegrate.replace(' ', '').split(',') if ',' in axesToIntegrate else [axesToIntegrate]
    integrate += axesToIntegrate
  hu = GetHisto(path, var, {'channel':ch, 'val':'num'}, group=['sample', 'pr', {'pr':mcSample}], integrate=integrate)
  hd = GetHisto(path, var, {'channel':ch, 'val':'den'}, group=['sample', 'pr', {'pr':mcSample}], integrate=integrate)
  du = GetHisto(path, var, {'channel':ch, 'val':'num'}, group=['sample', 'pr', {'pr':dataSample}], integrate=integrate)
  dd = GetHisto(path, var, {'channel':ch, 'val':'den'}, group=['sample', 'pr', {'pr':dataSample}], integrate=integrate)
  if var == 'pteta':
    ratio, do, up = GetSFfromCountsHisto(hu, hd, du, dd)
    cx = hu.axis('pt').edges()
    cy = hu.axis('abseta').edges()
    h2d = GetH2DfromXY([cx, cy], ratio, ytit='Y', xtit=['X1', 'X2'], hname=['pt','eta'])
    DrawEff2D(h2d,'pt', error=up, error2=do, xtit=xtit, ytit=ytit, tit=tit, outname=outname)
  else:
    DrawEff(hu, hd, du, dd, title=tit, xtit=xtit, doFill=None, mcColor='r', dataColor='k', outname=outname)


path = 'cafea/data/5TeV/triggerSFs.pkl.gz'
outpath = '/nfs/fanae/user/juanr/www/public/tt5TeV/ljets/SFs/'
DrawTriggerSFs(path, 'm', 'pteta', dataSample='HighEGJet', tit='Muon trigger scale factors', outname=outpath+'/MuonTrigSFs.png')
DrawTriggerSFs(path, 'e', 'pteta', dataSample='SingleMuon', tit='Electron trigger scale factors', outname=outpath+'/ElecTrigSFs.png')
DrawTriggerSFs(path, var='eta', ch='m', dataSample='HighEGJet', mcSample='tt', xtit='|$\eta$|', tit='Muon trigger efficiencies', outname=outpath+'/MuonTrigEff_eta.png')
DrawTriggerSFs(path, var='pt',  ch='m', dataSample='HighEGJet', mcSample='tt', xtit='p$_{T}$ (GeV)', tit='Muon trigger efficiencies', outname=outpath+'/MuonTrigEff_pt.png')
DrawTriggerSFs(path, var='eta', ch='e', dataSample='SingleMuon', mcSample='tt', xtit='|$\eta$|', tit='Electron trigger efficiencies', outname=outpath+'/ElecTrigEff_eta.png')
DrawTriggerSFs(path, var='pt', ch='e', dataSample='SingleMuon', mcSample='tt', xtit='p$_{T}$ (GeV)', tit='Electron trigger efficiencies', outname=outpath+'/ElecTrigEff_pt.png')
