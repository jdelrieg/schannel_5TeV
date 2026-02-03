# python ljets/test.py --path histos/2jun2022_btagM_jet25/
import copy
from config import *
from cafea.modules.fileReader import *

import statsmodels.api as sm
import argparse

if var is None: var = 'absu0eta'

plt = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi, var=var)
name = GetFileNameFromPath(path)


#hdampup,hdampdo = GetModSystHistos(path, 'TT_hdamp', 'hdamp', var=var)
#tuneup , tunedo = GetModSystHistos(path, 'TT_UE', 'UE', var=var)
#plt.AddExtraBkgHist([hdampup, hdampdo, tuneup, tunedo], add=True)

tuneup_top , tunedo_top = GetModSystHistos(path, 'tchannel_UE',    'UE', var=var)
tuneup_tbar , tunedo_tbar = GetModSystHistos(path, 'tbarchannel_UE',    'UE', var=var)
tuneup=tuneup_top+tuneup_tbar
tunedo=tunedo_top+tunedo_tbar

hdampup_top , hdampdo_top = GetModSystHistos(path, 'tchannel_hdamp',    'hdamp', var=var)
hdampup_tbar , hdampdo_tbar = GetModSystHistos(path, 'tbarchannel_hdamp',    'hdamp', var=var)
hdampup=hdampup_top+hdampup_tbar
hdampdo=hdampdo_top+hdampdo_tbar

plt.AddExtraBkgHist([hdampup,hdampdo,tuneup,tunedo],add=True)

'''
qcdpath=path+'shaped_QCD/'
QCDnom  = GetHisto(path+ 'QCD.pkl.gz', var, group=None)

QCD_shape_up= GetModSystHisto(qcdpath, 'QCD_splitshape_scaled',    'QCD_shapeUp', var=var,prname='QCD')
QCD_shape_do= GetModSystHisto(qcdpath, 'QCD_splitshape_scaled',    'QCD_shapeDown', var=var,prname='QCD')


QCD_shape_up_ll= GetModSystHisto(qcdpath, 'QCD_splitshape_scaled',    'QCD_shapeUp', var=var[:20]+'_ll'+'_cut',prname='QCD')
QCD_shape_up_lh=GetModSystHisto(qcdpath, 'QCD_splitshape_scaled',    'QCD_shapeUp', var=var[:20]+'_lh'+'_cut',prname='QCD')
QCD_shape_up_hl=GetModSystHisto(qcdpath, 'QCD_splitshape_scaled',    'QCD_shapeUp', var=var[:20]+'_hl'+'_cut',prname='QCD')
QCD_shape_up_hh=GetModSystHisto(qcdpath, 'QCD_splitshape_scaled',    'QCD_shapeUp', var=var[:20]+'_hh'+'_cut',prname='QCD')

QCD_shape_up.values()[('QCD',ch[0],level,'QCD_shapeUp')][()]=QCD_shape_up_ll.values()[('QCD',ch[0],level,'QCD_shapeUp')][()]+QCD_shape_up_lh.values()[('QCD',ch[0],level,'QCD_shapeUp')][()]+QCD_shape_up_hl.values()[('QCD',ch[0],level,'QCD_shapeUp')][()]+QCD_shape_up_hh.values()[('QCD',ch[0],level,'QCD_shapeUp')][()]

QCD_shape_do.values()[('QCD',ch[0],level,'QCD_shapeDown')][()]=QCDnom.values()[('QCD',ch[0],level,'norm')][()]-(QCD_shape_up.values()[('QCD',ch[0],level,'QCD_shapeUp')][()]-QCDnom.values()[('QCD',ch[0],level,'norm')][()])


plt.AddExtraBkgHist([QCD_shape_up,QCD_shape_do],add=True)
'''

plt.SetOutpath(path + '/combineFiles/')
plt.SetVerbose(verbose)
if not doData:
  plt.SetDataName('Asimov')
RebinVar(plt, var)

def SaveFile(level, channel):
  categories = {'level':level, 'channel':channel}
  lev  = categories['level'  ] if not isinstance(categories['level'  ], list) else categories['level'  ][0];
  chan = categories['channel'] if not isinstance(categories['channel'], list) else categories['channel'][0];
  outname = '%s_%s_%s'%(var, chan, lev)
  plt.SaveCombine(var, outname, categories=categories)

if __name__=="__main__":
  if not isinstance(ch, list): ch = [ch]
  if not isinstance(level, list): level = [level]
  for l in level:
    for c in ch:
      SaveFile(l, c)
