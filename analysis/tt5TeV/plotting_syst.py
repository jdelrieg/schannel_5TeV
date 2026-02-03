import gzip
import pickle
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
from coffea import hist, processor

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



# Load your data
with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/split_charge_goodJECs/tchannel_analysis.pkl.gz', 'rb') as file:
    loaded_data0 = pickle.load(file)
  
hep.style.use("CMS")

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/split_charge_goodJECs/tchannel_hdampUp.pkl.gz', 'rb') as file:
    loaded_hdampup = pickle.load(file)
with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/split_charge_goodJECs/tchannel_hdampDown.pkl.gz', 'rb') as file:
    loaded_hdampdown = pickle.load(file)

#with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/histos5TeV/tchannel_hdampDown_part1.pkl.gz', 'rb') as file:
#    loaded_hdampdown_extra = pickle.load(file)

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/split_charge_goodJECs/tchannel_UEUp.pkl.gz', 'rb') as file:
    loaded_ueup = pickle.load(file)
with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/split_charge_goodJECs/tchannel_UEDown.pkl.gz', 'rb') as file:
    loaded_uedown = pickle.load(file)



with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/split_charge_goodJECs/tbarchannel_analysis.pkl.gz', 'rb') as file:
    loaded_data0_bar = pickle.load(file)
  
hep.style.use("CMS")

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/split_charge_goodJECs/tbarchannel_hdampUp.pkl.gz', 'rb') as file:
    loaded_hdampup_bar = pickle.load(file)
with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/split_charge_goodJECs/tbarchannel_hdampDown.pkl.gz', 'rb') as file:
    loaded_hdampdown_bar = pickle.load(file)
with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/split_charge_goodJECs/tbarchannel_UEUp.pkl.gz', 'rb') as file:
    loaded_ueup_bar = pickle.load(file)
with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/split_charge_goodJECs/tbarchannel_UEDown.pkl.gz', 'rb') as file:
    loaded_uedown_bar = pickle.load(file)
    
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Description of your script')
	parser.add_argument('--uncs', action='store_true', help='Description of option')
	parser.add_argument('--doLowess', action='store_true', help='Description of option')
	parser.add_argument('--normerr', action='store_true', help='Description of option')
	parser.add_argument('--var', type=str, default='MVAscore_relaxed_b10', help='Description of option')
	parser.add_argument('--ch', type=str, default='e_minus', help='Description of option')
	parser.add_argument('--level', type=str, default='2j1b', help='Description of option')
	parser.add_argument('--nuisance', type=str, default='hdamp', help='Description of option')
	
	args = parser.parse_args()
	uncs=args.uncs
	doLowess=args.doLowess
	normerr=args.normerr
	var=args.var
	ch=args.ch
	level=args.level
	nuisance=args.nuisance	




        
	hnorm0,hnormerr0=loaded_data0[var].values(sumw2=True)[('tchannel', ch, level, 'norm')]
	h1=hnorm0
	h1error=hnormerr0
	h1=h1*302
	
	
#	h1=np.array([sum(h1[:4]),sum(h1[4:6]),sum(h1[6:8]),sum(h1[8:10]),sum(h1[10:12]),sum(h1[12:14]),sum(h1[14:16]),sum(h1[16:])]) if var=='MVAscore_relaxed_b10' else np.array(h1[:])#np.array([sum(h1[:5]),h1[5],h1[6],h1[7],h1[8],h1[9],sum(h1[10:])]) 
#	h1error=np.array([sum(h1error[:4]),sum(h1error[4:6]),sum(h1error[6:8]),sum(h1error[8:10]),sum(h1error[10:12]),sum(h1error[12:14]),sum(h1error[14:16]),sum(h1error[16:])]) if var=='MVAscore_relaxed_b10' else np.array(h1error[:])#np.array([sum(h1error[:5]),h1error[5],h1error[6],h1error[7],h1error[8],h1error[9],sum(h1error[10:])]) 
	
	hnorm0_bar,hnormerr0_bar=loaded_data0_bar[var].values(sumw2=True)[('tbarchannel', ch, level, 'norm')]
	h1_bar=hnorm0_bar
	h1error_bar=hnormerr0_bar
	h1_bar=h1_bar*302
	
	
	h1=np.array([sum(h1[:4])+sum(h1_bar[:4]),sum(h1[4:6])+sum(h1_bar[4:6]),sum(h1[6:8])+sum(h1_bar[6:8]),sum(h1[8:10])+sum(h1_bar[8:10]),sum(h1[10:12])+sum(h1_bar[10:12]),sum(h1[12:14])+sum(h1_bar[12:14]),sum(h1[14:16])+sum(h1_bar[14:16]),sum(h1[16:])+sum(h1_bar[16:])]) if var=='MVAscore_relaxed_b10' else np.array(h1[:]+h1_bar[:])#np.array([sum(h1[:5]),h1[5],h1[6],h1[7],h1[8],h1[9],sum(h1[10:])]) 
	h1error=np.array([sum(h1error[:4])+sum(h1error_bar[:4]),sum(h1error[4:6])+sum(h1error_bar[4:6]),sum(h1error[6:8])+sum(h1error_bar[6:8]),sum(h1error[8:10])+sum(h1error_bar[8:10]),sum(h1error[10:12])+sum(h1error_bar[10:12]),sum(h1error[12:14])+sum(h1error_bar[12:14]),sum(h1error[14:16])+sum(h1error_bar[14:16]),sum(h1error[16:])+sum(h1error_bar[16:])]) if var=='MVAscore_relaxed_b10' else np.array(h1error[:]+h1error_bar[:])#np.array([sum(h1error[:5]),h1error[5],h1error[6],h1error[7],h1error[8],h1error[9],sum(h1error[10:])]) 
#	h1=h1+h1_bar



	#For hdampUp    
	h2,h2error=loaded_hdampup[var].values(sumw2=True)[('tchannel_hdampUp',ch,level,'norm')] if nuisance=='hdamp' else loaded_ueup[var].values(sumw2=True)[('tchannel_UEUp',ch,level,'norm')]
	h2=h2*302

	h2_bar,h2error_bar=loaded_hdampup_bar[var].values(sumw2=True)[('tbarchannel_hdampUp',ch,level,'norm')] if nuisance=='hdamp' else loaded_ueup_bar[var].values(sumw2=True)[('tbarchannel_UEUp',ch,level,'norm')]
	h2_bar=h2_bar*302
	

	# up to here hdamp up


	#For hdampDown 
	h3,h3error=loaded_hdampdown[var].values(sumw2=True)[('tchannel_hdampDown',ch,level,'norm')] if nuisance=='hdamp' else loaded_uedown[var].values(sumw2=True)[('tchannel_UEDown',ch,level,'norm')]	
	h3=h3*302

#	h3_bis,h3error_bis=loaded_hdampdown_extra[var].values(sumw2=True)[('tchannel_hdampDown',ch,level,'norm')] if nuisance=='hdamp' else loaded_uedown[var].values(sumw2=True)[('tchannel_UEDown',ch,level,'norm')]	
#	h3_bis=h3_bis*302

	h3_bar,h3error_bar=loaded_hdampdown_bar[var].values(sumw2=True)[('tbarchannel_hdampDown',ch,level,'norm')] if nuisance=='hdamp' else loaded_uedown_bar[var].values(sumw2=True)[('tbarchannel_UEDown',ch,level,'norm')]	
	h3_bar=h3_bar*302

	#Up to here hdampdown
#	h2=h2+h2_bar
#	h3=h3+h3_bar
	h2=np.array([sum(h2[:4])+sum(h2_bar[:4]),sum(h2[4:6])+sum(h2_bar[4:6]),sum(h2[6:8])+sum(h2_bar[6:8]),sum(h2[8:10])+sum(h2_bar[8:10]),sum(h2[10:12])+sum(h2_bar[10:12]),sum(h2[12:14])+sum(h2_bar[12:14]),sum(h2[14:16])+sum(h2_bar[14:16]),sum(h2[16:])+sum(h2_bar[16:])]) if var=='MVAscore_relaxed_b10' else np.array(h2[:]+h2_bar[:])#np.array([sum(h2[:5]),h2[5],h2[6],h2[7],h2[8],h2[9],sum(h2[10:])]) 
	h2error=np.array([sum(h2error[:4])+sum(h2error_bar[:4]),sum(h2error[4:6])+sum(h2error_bar[4:6]),sum(h2error[6:8])+sum(h2error_bar[6:8]),sum(h2error[8:10])+sum(h2error_bar[8:10]),sum(h2error[10:12])+sum(h2error_bar[10:12]),sum(h2error[12:14])+sum(h2error_bar[12:14]),sum(h2error[14:16])+sum(h2error_bar[14:16]),sum(h2error[16:])+sum(h2error_bar[16:])]) if var=='MVAscore_relaxed_b10' else np.array(h2error[:]+h2error_bar[:])#np.array([sum(h2error[:5]),h2error[5],h2error[6],h2error[7],h2error[8],h2error[9],sum(h2error[10:])]) 
	h3=np.array([sum(h3[:4])+sum(h3_bar[:4]),sum(h3[4:6])+sum(h3_bar[4:6]),sum(h3[6:8])+sum(h3_bar[6:8]),sum(h3[8:10])+sum(h3_bar[8:10]),sum(h3[10:12])+sum(h3_bar[10:12]),sum(h3[12:14])+sum(h3_bar[12:14]),sum(h3[14:16])+sum(h3_bar[14:16]),sum(h3[16:])+sum(h3_bar[16:])]) if var=='MVAscore_relaxed_b10' else np.array(h3[:]+h3_bar[:])#np.array([sum(h3[:5]),h3[5],h3[6],h3[7],h3[8],h3[9],sum(h3[10:])]) 
	h3error=np.array([sum(h3error[:4])+sum(h3error_bar[:4]),sum(h3error[4:6])+sum(h3error_bar[4:6]),sum(h3error[6:8])+sum(h3error_bar[6:8]),sum(h3error[8:10])+sum(h3error_bar[8:10]),sum(h3error[10:12])+sum(h3error_bar[10:12]),sum(h3error[12:14])+sum(h3error_bar[12:14]),sum(h3error[14:16])+sum(h3error_bar[14:16]),sum(h3error[16:])+sum(h3error_bar[16:])]) if var=='MVAscore_relaxed_b10' else np.array(h3error[:]+h3error_bar[:])#np.array([sum(h3error[:5]),h3error[5],h3error[6],h3error[7],h3error[8],h3error[9],sum(h3error[10:])]) 

#	h3=np.array([sum(h3[:4])+sum(h3_bar[:4])+sum(h3_bis[:4]),sum(h3[4:6])+sum(h3_bar[4:6])+sum(h3_bis[4:6]),sum(h3[6:8])+sum(h3_bar[6:8])+sum(h3_bis[6:8]),sum(h3[8:10])+sum(h3_bar[8:10])+sum(h3_bis[8:10]),sum(h3[10:12])+sum(h3_bar[10:12])+sum(h3_bis[10:12]),sum(h3[12:14])+sum(h3_bar[12:14])+sum(h3_bis[12:14]),sum(h3[14:16])+sum(h3_bar[14:16])+sum(h3_bis[14:16]),sum(h3[16:])+sum(h3_bar[16:])+sum(h3_bis[16:])]) if var=='MVAscore_relaxed_b10' else np.array(h3[:]+h3_bar[:])#np.array([sum(h3[:5]),h3[5],h3[6],h3[7],h3[8],h3[9],sum(h3[10:])]) 
#	h3error=np.array([sum(h3error[:4])+sum(h3error_bar[:4]),sum(h3error[4:6])+sum(h3error_bar[4:6]),sum(h3error[6:8])+sum(h3error_bar[6:8]),sum(h3error[8:10])+sum(h3error_bar[8:10]),sum(h3error[10:12])+sum(h3error_bar[10:12]),sum(h3error[12:14])+sum(h3error_bar[12:14]),sum(h3error[14:16])+sum(h3error_bar[14:16]),sum(h3error[16:])+sum(h3error_bar[16:])]) if var=='MVAscore_relaxed_b10' else np.array(h3error[:]+h3error_bar[:])#np.array([sum(h3error[:5]),h3error[5],h3error[6],h3error[7],h3error[8],h3error[9],sum(h3error[10:])]) 

	

	
	bins = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5,6.5,7.5]) if var=='MVAscore_relaxed_b10' else np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5,6.5,7.5,8.5,9.5]) 
	bar_width = bins[1] - bins[0]
	if doLowess:
		
		xvalLowess = bins
		yvalUpLowess = h2
		yvalDnLowess = h3
		nominalLowess = h1
		
		ratio_diff = (yvalUpLowess/nominalLowess - yvalDnLowess/nominalLowess)/2.
		diff_smooth = sm.nonparametric.lowess(ratio_diff, xvalLowess)[:,1]

		# Chi2 part
		up_scale = get_smoothed_scale_factor( diff_smooth, nominalLowess, yvalUpLowess, np.var(nominalLowess), np.var(yvalUpLowess))
		down_scale = get_smoothed_scale_factor( diff_smooth, nominalLowess, yvalDnLowess, np.var(nominalLowess), np.var(yvalDnLowess))

		up_ratio = (1 + up_scale*diff_smooth)
		down_ratio = ( 1 + down_scale*diff_smooth)
		
		h2 = nominalLowess * np.nan_to_num(up_ratio, nan=1.0)
		h3 = nominalLowess * np.nan_to_num(down_ratio, nan=1.0)	


	# Plotting the ratio plot
	fig, (ax0, ax1) = plt.subplots(nrows=2, gridspec_kw={'height_ratios': [3, 1]}, sharex=True, figsize=(8, 8))





	hatch_pattern_up = '//'
	hatch_pattern_down = '\\'
	hatch_pattern_norm = '/\/\/'

	# Plot the upper values as bars
	if uncs==True:
		h2_statup=h2/302+np.sqrt(h2error)
		h2_statdown=h2/302-np.sqrt(h2error)
		h2_statup=h2_statup*302
		h2_statdown=h2_statdown*302
		
	#	h2_statup=np.array([sum(h2_statup[:6]),sum(h2_statup[6:8]),sum(h2_statup[8:10]),sum(h2_statup[10:12]),sum(h2_statup[12:14]),sum(h2_statup[14:])])
	#	h2_statdown=np.array([sum(h2_statdown[:6]),sum(h2_statdown[6:8]),sum(h2_statdown[8:10]),sum(h2_statdown[10:12]),sum(h2_statdown[12:14]),sum(h2_statdown[14:])])
		
		h2_statup_ratio=h2_statup/h1
		h2_statdown_ratio=h2_statdown/h1
		
		h3_statup=h3/302+np.sqrt(h3error)
		h3_statdown=h3/302-np.sqrt(h3error)
		h3_statup=h3_statup*302
		h3_statdown=h3_statdown*302
		
	#	h3_statup=np.array([sum(h3_statup[:6]),sum(h3_statup[6:8]),sum(h3_statup[8:10]),sum(h3_statup[10:12]),sum(h3_statup[12:14]),sum(h3_statup[14:])])
	#	h3_statdown=np.array([sum(h3_statdown[:6]),sum(h3_statdown[6:8]),sum(h3_statdown[8:10]),sum(h3_statdown[10:12]),sum(h3_statdown[12:14]),sum(h3_statdown[14:])])
		
		h3_statup_ratio=h3_statup/h1
		h3_statdown_ratio=h3_statdown/h1
		
		ax0.bar(bins, h2_statup - h2_statdown, bottom=h2_statdown, width=bar_width, color='none', hatch=hatch_pattern_up, edgecolor='red',linewidth=0.2, label='Up Stat unc')
		ax0.bar(bins, h3_statup - h3_statdown, bottom=h3_statdown, width=bar_width, color='none', hatch=hatch_pattern_down, edgecolor='blue',linewidth=0.2, label='Down Stat unc')
		ax1.bar(bins, h2_statup_ratio - h2_statdown_ratio, bottom=h2_statdown_ratio, width=bar_width, color='none', hatch=hatch_pattern_up, edgecolor='red',linewidth=0.5, label='Up Stat unc')
		ax1.bar(bins, h3_statup_ratio - h3_statdown_ratio, bottom=h3_statdown_ratio, width=bar_width, color='none', hatch=hatch_pattern_down, edgecolor='blue',linewidth=0.5, label='Down Stat unc')

	if normerr: 
		h1_statup=h1/302+np.sqrt(h1error)
		h1_statdown=h1/302-np.sqrt(h1error)
		h1_statup=h1_statup*302
		h1_statdown=h1_statdown*302	
		ax0.bar(bins, h1_statup - h1_statdown, bottom=h1_statdown, width=bar_width, color='none', hatch=hatch_pattern_norm, edgecolor='gray',linewidth=0.2, label='Stat unc')
		statup_ratio=h1_statup/h1
		statdown_ratio=h1_statdown/h1
		ax1.bar(bins, statup_ratio - statdown_ratio, bottom=statdown_ratio, width=bar_width, color='none', hatch=hatch_pattern_norm, alpha=0.5,edgecolor='gray',linewidth=0, label='Stat unc')

	hep.histplot(h1, ax=ax0, stack=False, histtype="step", color=['black'], label=['norm'])
	hep.histplot(h2, ax=ax0, stack=False, histtype="step", color=['red'], label=['Up'])
	hep.histplot(h3, ax=ax0, stack=False, histtype="step", color=['blue'], label=['Down'])

	ax0.legend(fontsize=10,ncol=3)

	ratioup=h2/h1
	hep.histplot(ratioup, ax=ax1, color=['red'], histtype='step')
	#if uncs==True: ax1.bar(bins, h2_statup_ratio - h2_statdown_ratio, bottom=h2_statdown_ratio, width=bar_width, color='none', hatch=hatch_pattern_up, edgecolor='red',linewidth=0.2, label='Up Stat unc')

	ratiodown=h3/h1
	hep.histplot(ratiodown, ax=ax1, color=['blue'], histtype='step')
	#if uncs==True: ax1.bar(bins, h3_statup_ratio - h3_statdown_ratio, bottom=h3_statdown_ratio, width=bar_width, color='none', hatch=hatch_pattern_down, edgecolor='blue',linewidth=0.2, label='Down Stat unc')


	ax1.axhline(y=1, color='gray', linestyle='--')

	ax1.set_ylim(0.85, 1.15)

	# Set labels and title
	ax0.set_ylabel('Events')
	ax1.set_xlabel(var)
	ax1.set_ylabel('Ratio')
	ax0.text(0.05,0.9,ch+level+', '+nuisance,fontsize=15,transform=ax0.transAxes)
	#plt.savefig("ratio_scale_3j1b{}.png".format(num_pdfs-1))
	plt.savefig("lowess/{}_{}_{}{}.png".format(nuisance,'no_lowess' if doLowess==False else 'lowess', ch+level,'_uncs' if uncs else None))
	exit()



