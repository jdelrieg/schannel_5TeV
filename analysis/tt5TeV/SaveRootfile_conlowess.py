# python ljets/test.py --path histos/2jun2022_btagM_jet25/
import copy
from config import *
from cafea.modules.fileReader import *

import statsmodels.api as sm
import argparse

def get_smoothed_scale_factor( smoothed_ratio_diff, nom_hist, var_hist, nom_var, var_var):
    '''determine an overall systematic template smoothing which minimizes the chi-squared between
        the smoothed systematic and the unsmoothed template it is derived from.
        Taken from:
             https://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2018_077_v4.pdf'''

    total_var = np.sqrt( nom_var + var_var)
    unsmoothed_diff = var_hist - nom_hist
    ratio_factor = smoothed_ratio_diff * nom_hist
    numerator = np.sum( ratio_factor * unsmoothed_diff / total_var )
    denominator = np.sum( (ratio_factor / np.sqrt(total_var))**2 )
    return numerator/denominator


if var is None: var = 'absu0eta'

plt = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi, var=var)
name = GetFileNameFromPath(path)

htop  = GetHisto(path+ 'tchannel_analysis' +   '.pkl.gz', var, group=None)
htbar=GetHisto(path+ 'tbarchannel_analysis' +   '.pkl.gz', var, group=None) #These two lines load the nominal



tuneup_top , tunedo_top = GetModSystHistos(path, 'tchannel_UE',    'UE', var=var)
tuneup_tbar , tunedo_tbar = GetModSystHistos(path, 'tbarchannel_UE',    'UE', var=var)
#tuneup=tuneup_top+tuneup_tbar
#tunedo=tunedo_top+tunedo_tbar #These 4 lines load the classical UE histos



tunehistoup_top=copy.deepcopy(tuneup_top)
tunehistoup_tbar=copy.deepcopy(tuneup_tbar)
data_dict_up_top=tunehistoup_top.values()
data_dict_up_tbar=tunehistoup_tbar.values()
tunehistodo_top=copy.deepcopy(tunedo_top)
tunehistodo_tbar=copy.deepcopy(tunedo_tbar)
data_dict_do_top=tunehistodo_top.values()
data_dict_do_tbar=tunehistodo_tbar.values() #We make copies and build an object type that allow us making changes and saving them


#Lowess starting for UE
#xvalLowess = np.array([0.15, 0.25, 0.35, 0.45, 0.55, 0.65,0.75,0.85]) if var=='MVAscore_relaxed_b10' else np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75,3.25,3.75,4.25,4.75]) 
xvalLowess = np.array([0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425,0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875,0.925, 0.975]) if var=='MVAscore_relaxed_b10' else np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75,3.25,3.75,4.25,4.75]) 
#xvalLowess = np.array([ 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875,0.925, 0.975]) if var=='MVAscore_relaxed_b10' else np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75,3.25,3.75,4.25,4.75])

yvalUpLowess_top = data_dict_up_top[('tchan',ch[0],level,'UEUp')]
yvalDnLowess_top = data_dict_do_top[('tchan',ch[0],level,'UEDown')]
nominalLowess_top = htop.values()[('tchannel', ch[0], level, 'norm')]

yvalUpLowess_tbar = data_dict_up_tbar[('tbarchan',ch[0],level,'UEUp')]
yvalDnLowess_tbar = data_dict_do_tbar[('tbarchan',ch[0],level,'UEDown')]
nominalLowess_tbar=htbar.values()[('tbarchannel', ch[0], level, 'norm')]


ratio_diff_top = (yvalUpLowess_top/nominalLowess_top - yvalDnLowess_top/nominalLowess_top)/2.
ratio_diff_top = np.nan_to_num(ratio_diff_top, nan=0)
diff_smooth_top = sm.nonparametric.lowess(ratio_diff_top, xvalLowess)[:,1]

ratio_diff_tbar = (yvalUpLowess_tbar/nominalLowess_tbar - yvalDnLowess_tbar/nominalLowess_tbar)/2.
ratio_diff_tbar = np.nan_to_num(ratio_diff_tbar, nan=0)
diff_smooth_tbar = sm.nonparametric.lowess(ratio_diff_tbar, xvalLowess)[:,1]


# Chi2 part
up_scale_top = get_smoothed_scale_factor( diff_smooth_top, nominalLowess_top, yvalUpLowess_top, np.var(nominalLowess_top), np.var(yvalUpLowess_top))
down_scale_top = get_smoothed_scale_factor( diff_smooth_top, nominalLowess_top, yvalDnLowess_top, np.var(nominalLowess_top), np.var(yvalDnLowess_top))

up_ratio_top = (1 + up_scale_top*diff_smooth_top)
down_ratio_top = ( 1 + down_scale_top*diff_smooth_top)


up_scale_tbar = get_smoothed_scale_factor( diff_smooth_tbar, nominalLowess_tbar, yvalUpLowess_tbar, np.var(nominalLowess_tbar), np.var(yvalUpLowess_tbar))
down_scale_tbar = get_smoothed_scale_factor( diff_smooth_tbar, nominalLowess_tbar, yvalDnLowess_tbar, np.var(nominalLowess_tbar), np.var(yvalDnLowess_tbar))

up_ratio_tbar = (1 + up_scale_tbar*diff_smooth_tbar)
down_ratio_tbar = ( 1 + down_scale_tbar*diff_smooth_tbar)


tunehistoup_top.values()[('tchan',ch[0],level,'UEUp')][()]= nominalLowess_top * np.nan_to_num(up_ratio_top, nan=1.0)
tunehistodo_top.values()[('tchan',ch[0],level,'UEDown')][()] = nominalLowess_top * np.nan_to_num(down_ratio_top, nan=1.0)

tunehistoup_tbar.values()[('tbarchan',ch[0],level,'UEUp')][()]= nominalLowess_tbar * np.nan_to_num(up_ratio_tbar, nan=1.0)
tunehistodo_tbar.values()[('tbarchan',ch[0],level,'UEDown')][()] = nominalLowess_tbar * np.nan_to_num(down_ratio_tbar, nan=1.0)	#These two lines assign the lowess-returned values (right hand side of eq) to the copies of the histos we made at the beginning (left hand side of eq)
#Lowess up to here for UE

del(xvalLowess,yvalUpLowess_top,yvalDnLowess_top,nominalLowess_top,ratio_diff_top,diff_smooth_top,up_scale_top,down_scale_top,up_ratio_top,down_ratio_top,yvalUpLowess_tbar,yvalDnLowess_tbar,nominalLowess_tbar,ratio_diff_tbar,diff_smooth_tbar,up_scale_tbar,down_scale_tbar,up_ratio_tbar,down_ratio_tbar)

#And now we repeat for hdamp
hdampup_top , hdampdo_top = GetModSystHistos(path, 'tchannel_hdamp',    'hdamp', var=var)
hdampup_tbar , hdampdo_tbar = GetModSystHistos(path, 'tbarchannel_hdamp',    'hdamp', var=var)
#hdampup=hdampup_top+hdampup_tbar
#hdampdo=hdampdo_top+hdampdo_tbar

hdamphistoup_top=copy.deepcopy(hdampup_top)
data_dict_up_hdamp_top=hdamphistoup_top.values()
hdamphistodo_top=copy.deepcopy(hdampdo_top)
data_dict_do_hdamp_top=hdamphistodo_top.values() 

hdamphistoup_tbar=copy.deepcopy(hdampup_tbar)
data_dict_up_hdamp_tbar=hdamphistoup_tbar.values()
hdamphistodo_tbar=copy.deepcopy(hdampdo_tbar)
data_dict_do_hdamp_tbar=hdamphistodo_tbar.values() #We make copies and build an object type that allow us making changes and saving them


#Lowess starting for UE
#xvalLowess = np.array([0.15, 0.25, 0.35, 0.45, 0.55, 0.65,0.75,0.85]) if var=='MVAscore_relaxed_b10' else np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75,3.25,3.75,4.25,4.75]) 
xvalLowess = np.array([0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425,0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875,0.925, 0.975]) if var=='MVAscore_relaxed_b10' else np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75,3.25,3.75,4.25,4.75]) 
#xvalLowess = np.array([ 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875,0.925, 0.975]) if var=='MVAscore_relaxed_b10' else np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75,3.25,3.75,4.25,4.75])

yvalUpLowess_top = data_dict_up_hdamp_top[('tchan',ch[0],level,'hdampUp')]
yvalDnLowess_top = data_dict_do_hdamp_top[('tchan',ch[0],level,'hdampDown')]
nominalLowess_top = htop.values()[('tchannel', ch[0], level, 'norm')]

yvalUpLowess_tbar = data_dict_up_hdamp_tbar[('tbarchan',ch[0],level,'hdampUp')]
yvalDnLowess_tbar = data_dict_do_hdamp_tbar[('tbarchan',ch[0],level,'hdampDown')]
nominalLowess_tbar=htbar.values()[('tbarchannel', ch[0], level, 'norm')]


ratio_diff_top = (yvalUpLowess_top/nominalLowess_top - yvalDnLowess_top/nominalLowess_top)/2.
ratio_diff_top = np.nan_to_num(ratio_diff_top, nan=0)
diff_smooth_top = sm.nonparametric.lowess(ratio_diff_top, xvalLowess)[:,1]

ratio_diff_tbar = (yvalUpLowess_tbar/nominalLowess_tbar - yvalDnLowess_tbar/nominalLowess_tbar)/2.
ratio_diff_tbar = np.nan_to_num(ratio_diff_tbar, nan=0)
diff_smooth_tbar = sm.nonparametric.lowess(ratio_diff_tbar, xvalLowess)[:,1]


# Chi2 part
up_scale_top = get_smoothed_scale_factor( diff_smooth_top, nominalLowess_top, yvalUpLowess_top, np.var(nominalLowess_top), np.var(yvalUpLowess_top))
down_scale_top = get_smoothed_scale_factor( diff_smooth_top, nominalLowess_top, yvalDnLowess_top, np.var(nominalLowess_top), np.var(yvalDnLowess_top))
up_ratio_top = (1 + up_scale_top*diff_smooth_top)
down_ratio_top = ( 1 + down_scale_top*diff_smooth_top)

up_scale_tbar = get_smoothed_scale_factor( diff_smooth_tbar, nominalLowess_tbar, yvalUpLowess_tbar, np.var(nominalLowess_tbar), np.var(yvalUpLowess_tbar))
down_scale_tbar = get_smoothed_scale_factor( diff_smooth_tbar, nominalLowess_tbar, yvalDnLowess_tbar, np.var(nominalLowess_tbar), np.var(yvalDnLowess_tbar))
up_ratio_tbar = (1 + up_scale_tbar*diff_smooth_tbar)
down_ratio_tbar = ( 1 + down_scale_tbar*diff_smooth_tbar)


hdamphistoup_top.values()[('tchan',ch[0],level,'hdampUp')][()]= nominalLowess_top * np.nan_to_num(up_ratio_top, nan=1.0)
hdamphistodo_top.values()[('tchan',ch[0],level,'hdampDown')][()] = nominalLowess_top * np.nan_to_num(down_ratio_top, nan=1.0)

hdamphistoup_tbar.values()[('tbarchan',ch[0],level,'hdampUp')][()]= nominalLowess_tbar * np.nan_to_num(up_ratio_tbar, nan=1.0)
hdamphistodo_tbar.values()[('tbarchan',ch[0],level,'hdampDown')][()] = nominalLowess_tbar * np.nan_to_num(down_ratio_tbar, nan=1.0)	#These two lines assign the lowess-returned values (right hand side of eq) to the copies of the histos we made at the beginning (left hand side of eq)
#Lowess up to here for UE

del(xvalLowess,yvalUpLowess_top,yvalDnLowess_top,nominalLowess_top,ratio_diff_top,diff_smooth_top,up_scale_top,down_scale_top,up_ratio_top,down_ratio_top,yvalUpLowess_tbar,yvalDnLowess_tbar,nominalLowess_tbar,ratio_diff_tbar,diff_smooth_tbar,up_scale_tbar,down_scale_tbar,up_ratio_tbar,down_ratio_tbar)


plt.AddExtraBkgHist([hdamphistoup_top,hdamphistodo_top,hdamphistoup_tbar,hdamphistodo_tbar,tunehistoup_top,tunehistodo_top,tunehistoup_tbar,tunehistodo_tbar],add=True) #If lowess
#plt.AddExtraBkgHist([hdampup,hdampdo,tuneup,tunedo],add=True) #If no lowess





     #Second attempt of shapes
qcdpath=path+'shaped_QCD/'
QCDnom  = GetHisto(path+ 'QCD.pkl.gz', var, group=None)

QCD_shape_up= GetModSystHisto(qcdpath, 'QCD_splitshape_scaled',    'QCD_shapeUp', var=var,prname='QCD')
QCD_shape_do= GetModSystHisto(qcdpath, 'QCD_splitshape_scaled',    'QCD_shapeDown', var=var,prname='QCD')

QCD_shape_up_ll= GetModSystHisto(qcdpath, 'QCD_splitshape_scaled',    'QCD_shapeUp', var=var+'_ll',prname='QCD')
QCD_shape_up_lh=GetModSystHisto(qcdpath, 'QCD_splitshape_scaled',    'QCD_shapeUp', var=var+'_lh',prname='QCD')
QCD_shape_up_hl=GetModSystHisto(qcdpath, 'QCD_splitshape_scaled',    'QCD_shapeUp', var=var+'_hl',prname='QCD')
QCD_shape_up_hh=GetModSystHisto(qcdpath, 'QCD_splitshape_scaled',    'QCD_shapeUp', var=var+'_hh',prname='QCD')

QCD_shape_up.values()[('QCD',ch[0],level,'QCD_shapeUp')][()]=QCD_shape_up_ll.values()[('QCD',ch[0],level,'QCD_shapeUp')][()]+QCD_shape_up_lh.values()[('QCD',ch[0],level,'QCD_shapeUp')][()]+QCD_shape_up_hl.values()[('QCD',ch[0],level,'QCD_shapeUp')][()]+QCD_shape_up_hh.values()[('QCD',ch[0],level,'QCD_shapeUp')][()]

QCD_shape_do.values()[('QCD',ch[0],level,'QCD_shapeDown')][()]=np.maximum(QCDnom.values()[('QCD',ch[0],level,'norm')][()]-(QCD_shape_up.values()[('QCD',ch[0],level,'QCD_shapeUp')][()]-QCDnom.values()[('QCD',ch[0],level,'norm')][()]),0)
															   #Nopte that since we are simmetrizing, it may happen that there are negative values. We set them to 0 if that happens with the np.maximum(value,0) that we do in this line

plt.AddExtraBkgHist([QCD_shape_up,QCD_shape_do],add=True)




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

######################################################################################## DOWN HERE IT WAS THE FIRST ATTEMP FOR QCD SHAPES. NO LONGER USED BUT KEPT COMMENTED IN THE CASE IT IS NEEDED FOR THE FUTURE

'''       First attempt of shapes
##Adding QCD ##
qcdpath=path+'binnedQCD/'
QCDnom  = GetHisto(path+ 'QCD.pkl.gz', var, group=None)

QCD_ll_up= GetModSystHisto(qcdpath, 'QCD_ll',    'QCD_llUp', var=var,prname='QCD')
QCD_ll_do= GetModSystHisto(qcdpath, 'QCD_ll',    'QCD_llDown', var=var,prname='QCD')

QCD_ll_do.values()[('QCD',ch[0],level,'QCD_llDown')][()]=QCDnom.values()[('QCD',ch[0],level,'norm')][()]-(QCD_ll_up.values()[('QCD',ch[0],level,'QCD_llUp')][()]-QCDnom.values()[('QCD',ch[0],level,'norm')][()])

QCD_lh_up= GetModSystHisto(qcdpath, 'QCD_lh',    'QCD_lhUp', var=var,prname='QCD')
QCD_lh_do= GetModSystHisto(qcdpath, 'QCD_lh',    'QCD_lhDown', var=var,prname='QCD')
QCD_lh_do.values()[('QCD',ch[0],level,'QCD_lhDown')][()]=QCDnom.values()[('QCD',ch[0],level,'norm')][()]-(QCD_lh_up.values()[('QCD',ch[0],level,'QCD_lhUp')][()]-QCDnom.values()[('QCD',ch[0],level,'norm')][()])

QCD_hl_up= GetModSystHisto(qcdpath, 'QCD_hl',    'QCD_hlUp', var=var,prname='QCD')
QCD_hl_do= GetModSystHisto(qcdpath, 'QCD_hl',    'QCD_hlDown', var=var,prname='QCD')
QCD_hl_do.values()[('QCD',ch[0],level,'QCD_hlDown')][()]=QCDnom.values()[('QCD',ch[0],level,'norm')][()]-(QCD_hl_up.values()[('QCD',ch[0],level,'QCD_hlUp')][()]-QCDnom.values()[('QCD',ch[0],level,'norm')][()])

QCD_hh_up= GetModSystHisto(qcdpath, 'QCD_hh',    'QCD_hhUp', var=var,prname='QCD')
QCD_hh_do= GetModSystHisto(qcdpath, 'QCD_hh',    'QCD_hhDown', var=var,prname='QCD')
QCD_hh_do.values()[('QCD',ch[0],level,'QCD_hhDown')][()]=QCDnom.values()[('QCD',ch[0],level,'norm')][()]-(QCD_hh_up.values()[('QCD',ch[0],level,'QCD_hhUp')][()]-QCDnom.values()[('QCD',ch[0],level,'norm')][()])

plt.AddExtraBkgHist([QCD_ll_up,QCD_ll_do,QCD_lh_up,QCD_lh_do,QCD_hl_up,QCD_hl_do,QCD_hh_up,QCD_hh_do],add=True)
###### End adding QCD
'''
