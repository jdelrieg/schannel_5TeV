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
from cycler import cycler
from topcoffea.plotter.OutText import OutText
from config import *

# Select variable, channel and cuts
year=2000
var = 'counts'
ch=['m']
le='incl' #'g2jets'#dilep
sy='norm'
#outname='/mnt_pool/c3_users/user/andreatf/www/private/ttrun3/stack_2018/%s_%s_%s.png'%(var,ch[0],le)
outname = '/nfs/fanae/user/andreatf/www/public/tt5TeV/checkNanoVicor/%s_%s_%s_hoyhoy.png'%(var,ch[0],le)
path='kk/TTPS.pkl.gz'
#path1=','histos5TeV_splitJES_0b/W2JetsToLNu_part1.pkl.gz','histos5TeV_splitJES_0b/W3JetsToLNu_part1.pkl.gz','histos5TeV_splitJES_0b/W2JetsToLNu_part2.pkl.gz','histos5TeV_splitJES_0b/W3JetsToLNu_part2.pkl.gz','histos5TeV_splitJES_0b/W2JetsToLNu_part3.pkl.gz','histos5TeV_splitJES_0b/W3JetsToLNu_part3.pkl.gz']
path1='kk/TTPS_vicorhoyhoy.pkl.gz'
def getHist(path,isData=False,overlay=False,isFake=False):
	hists = {}
	listpath = path if isinstance(path, list) else [path]
	for path in listpath:
		with gzip.open(path) as fin:
			hin = pickle.load(fin)
			for k in hin.keys():
				if k in hists: hists[k]+=hin[k]
				else:               hists[k]=hin[k]
	h = hists[var]
	h = h.integrate('channel', ch)
	h = h.integrate('level', le)
	h = h.integrate('syst',sy)
	if overlay==False: h = h.sum('sample')
	if isData==False: 
		if year==2018: h.scale(1000*59.83)
		elif year==2017: h.scale(1000*41.48)
		elif year==2016: h.scale(1000*16.81)
		elif year==2022: h.scale(1000*1.2)
		else: h.scale(302)
	return(h)

h1=getHist(path, overlay=False)
h2=getHist(path1, overlay=False)

'''
h1.scale(1./h1.values()[()].sum())
h2.scale(1./h2.values()[()].sum())
'''
if var=='met':
  h2=h2.rebin('met', hist.Bin("newmet","MET (GeV)", 20, 0, 200))
  h1=h1.rebin('met', hist.Bin("newmet","MET (GeV)", 20, 0, 200))

fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
plt.subplots_adjust(hspace=.07)
#fig = plt.figure(figsize=(14,7))

dataOpts = {'linestyle':'none', 'marker':'.', 'markersize':10., 'color':'k', 'elinewidth':1}

# Plot and save figure to outname
#hist.plot1d(h_bkg,ax=ax,clear=False,stack=True)#,fill_opts={'color':'green'})#
hist.plot1d(h1,ax=ax,clear=False,stack=False,line_opts={'color':'red','linewidth':3})#,overlay='sample'
hist.plot1d(h2,ax=ax,clear=False,stack=False,line_opts={'color':'blue','linewidth':3})#,overlay='sample'
#hist.plot1d(h3,ax=ax,clear=False,stack=False,line_opts={'color':'orange','linewidth':3})#,overlay='sample'
#hist.plot1d(h4,ax=ax,clear=False,stack=False,line_opts={'color':'green','linewidth':3})#,overlay='sample'

plt.text(5., 100000, r"$\bf{CMS}$ Preliminary", fontsize='x-large', horizontalalignment='left', verticalalignment='bottom')#, transform=ax.transAxes)
plt.text(3., 100000, r"%1.0f %s (%s)"%(1200, 'pb', '13 TeV'), fontsize='x-large', horizontalalignment='right', verticalalignment='bottom')#, transform=ax.transAxes)


ax.legend(('Xuan','Vicor hoy'),fontsize='x-large')
#ax.legend(('PFJets','PuppiJets'),fontsize='x-large')


hist.plotratio(
    num=h1,
    denom=h2,
    ax=rax,
    error_opts=dataOpts,
    denom_fill_opts={},
    guide_opts={},
    unc='num'
)
plt.ylabel('Ratio',fontsize='x-large')
plt.ylim(0.5,1.5)
#plt.xlim(80,105)
fig.savefig(outname)
print('Output histogram saved in %s'%outname)


