import pickle
import gzip
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
from coffea import hist, processor


with gzip.open('/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/split_charge_goodJECs_mistag_comb_btagEff/QCD_shape/QCD_shapes.pkl.gz','rb') as file: #
	loaded_data=pickle.load(file)

with gzip.open('/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/QCD_forEVAN_AN/QCD_removing_ht_mtw.pkl.gz','rb') as file: #
	loaded_data1=pickle.load(file)
	
with gzip.open('/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/QCD_forEVAN_AN/QCD_removing_ht.pkl.gz','rb') as file: #
	loaded_data2=pickle.load(file)
	
hep.style.use("CMS")

var='MVAscore_relaxed_b10'
channels=['e_fake_minus','e_fake_plus','m_fake_minus','m_fake_plus']
levels=['2j1b','3j1b','3j2b']

h0=np.zeros(len(loaded_data[var].values()[('QCD', 'm_fake_minus', '2j1b', 'norm')]))
p0=np.zeros(len(loaded_data[var].values()[('QCD', 'm_fake_minus', '2j1b', 'norm')]))
p1=np.zeros(len(loaded_data[var].values()[('QCD', 'm_fake_minus', '2j1b', 'norm')]))
p2=np.zeros(len(loaded_data[var].values()[('QCD', 'm_fake_minus', '2j1b', 'norm')]))
h1=np.zeros(len(loaded_data1[var].values()[('QCD', 'm_fake_minus', '2j1b', 'norm')]))
h2=np.zeros(len(loaded_data1[var].values()[('QCD', 'm_fake_minus', '2j1b', 'norm')]))

for ch in channels:
	for level in levels:
		if ch in ['e_fake_minus','e_fake_plus']:
			if (ch=='e_fake_plus')&(level=='3j2b'):
				print('canla') 
			else:
				 h0+=loaded_data[var].values()[('QCD', ch, level, 'norm')]
				 p0+=loaded_data[var].values(sumw2=True)[('QCD', ch, level, 'norm')][1]
		else:
			h0+=loaded_data[var].values()[('QCD', ch, level, 'norm')]
			p0+=loaded_data[var].values(sumw2=True)[('QCD', ch, level, 'norm')][1] 

for ch in channels:
	for level in levels:
		if ch in ['e_fake_minus','e_fake_plus']:
			if ((ch=='e_fake_minus')&(level=='2j1b'))|((ch=='e_fake_minus')&(level=='3j2b'))|((ch=='e_fake_plus')&(level=='2j1b'))|((ch=='e_fake_plus')&(level=='3j1b'))|((ch=='e_fake_plus')&(level=='3j2b')):
				print('canla') 
			else:
				 h1+=loaded_data1[var].values()[('QCD', ch, level, 'norm')]
				 p1+=loaded_data1[var].values(sumw2=True)[('QCD', ch, level, 'norm')][1]
		else:
			h1+=loaded_data1[var].values()[('QCD', ch, level, 'norm')]
			p1+=loaded_data1[var].values(sumw2=True)[('QCD', ch, level, 'norm')][1] 


for ch in channels:
	for level in levels:
		if ch in ['e_fake_minus','e_fake_plus']:
			if (ch=='e_fake_minus')|(ch=='e_fake_plus'):
				print('canla') 
			else:
				 h2+=loaded_data2[var].values()[('QCD', ch, level, 'norm')]
				 p2+=loaded_data2[var].values(sumw2=True)[('QCD', ch, level, 'norm')][1]
		else:
			if ((ch=='m_fake_plus')&(level=='3j2b')):
				print('canal')
			else:
				 h2+=loaded_data2[var].values()[('QCD', ch, level, 'norm')]
				 p2+=loaded_data2[var].values(sumw2=True)[('QCD', ch, level, 'norm')][1] 

'''
for ch in channels:
	for level in levels:
		if ch in ['e_fake_minus','e_fake_plus']:
			if (ch=='e_fake_plus')&(level=='3j2b_'):
				print('canla') 
			else:																#TRICK IF PLOTTING E/M ONLY VARIABLES
				 var='ept'
				 h0+=loaded_data[var].values()[('QCD', ch, level, 'norm')]
		else:
			var='mpt'
			h0+=loaded_data[var].values()[('QCD', ch, level, 'norm')] 
'''
plt.figure(figsize=(12, 8))  # Define the figure size


unc0up=h0+np.sqrt(p0)
unc0do=h0-np.sqrt(p0)

unc1up=h1+np.sqrt(p1)
unc1do=h1-np.sqrt(p1)

unc2up=h2+np.sqrt(p2)
unc2do=h2-np.sqrt(p2)





if var=='MVAscore_relaxed_b10':
	h1=np.array([sum(h1[:4]),sum(h1[4:6]),sum(h1[6:8]),sum(h1[8:10]),sum(h1[10:12]),sum(h1[12:14]),sum(h1[14:])])#]),sum(h1[16:])]) 
	h0=np.array([sum(h0[:4]),sum(h0[4:6]),sum(h0[6:8]),sum(h0[8:10]),sum(h0[10:12]),sum(h0[12:14]),sum(h0[14:])])#,sum(h0[16:])]) 
	h2=np.array([sum(h2[:4]),sum(h2[4:6]),sum(h2[6:8]),sum(h2[8:10]),sum(h2[10:12]),sum(h2[12:14]),sum(h2[14:])])#,sum(h2[16:])]) 
	
	unc0up=np.array([sum(unc0up[:4]),sum(unc0up[4:6]),sum(unc0up[6:8]),sum(unc0up[8:10]),sum(unc0up[10:12]),sum(unc0up[12:14]),sum(unc0up[14:])])#,sum(unc0up[16:])]) 
	unc0do=np.array([sum(unc0do[:4]),sum(unc0do[4:6]),sum(unc0do[6:8]),sum(unc0do[8:10]),sum(unc0do[10:12]),sum(unc0do[12:14]),sum(unc0do[14:])])#,sum(unc0do[16:])]) 
	unc1up=np.array([sum(unc1up[:4]),sum(unc1up[4:6]),sum(unc1up[6:8]),sum(unc1up[8:10]),sum(unc1up[10:12]),sum(unc1up[12:14]),sum(unc1up[14:])])#,sum(unc1up[16:])]) 
	unc1do=np.array([sum(unc1do[:4]),sum(unc1do[4:6]),sum(unc1do[6:8]),sum(unc1do[8:10]),sum(unc1do[10:12]),sum(unc1do[12:14]),sum(unc1do[14:])])#,sum(unc1do[16:])]) 
	unc2up=np.array([sum(unc2up[:4]),sum(unc2up[4:6]),sum(unc2up[6:8]),sum(unc2up[8:10]),sum(unc2up[10:12]),sum(unc2up[12:14]),sum(unc2up[14:])])#,sum(unc2up[16:])]) 
	unc2do=np.array([sum(unc2do[:4]),sum(unc2do[4:6]),sum(unc2do[6:8]),sum(unc2do[8:10]),sum(unc2do[10:12]),sum(unc2do[12:14]),sum(unc2do[14:])])#,sum(unc2do[16:])]) 
	
variacion2=unc2up/h2-np.ones(len(h2))
h2=h2/sum(h2)
unc2up=h2+h2*variacion2
unc2do=h2-h2*variacion2

variacion1=unc1up/h1-np.ones(len(h1))
h1=h1/sum(h1)
unc1up=h1+h1*variacion1
unc1do=h1-h1*variacion1

variacion0=unc0up/h0-np.ones(len(h0))
h0=h0/sum(h0)
unc0up=h0+h0*variacion0
unc0do=h0-h0*variacion0

# Plot the histogram of h0 with some customization
bin_edges = np.array([0, 1, 2, 3, 4, 5,6,7]) 
bin_edges_end=bin_edges+np.ones(len(bin_edges))

#plt.hlines(h0,bin_edges[:-1],bin_edges_end[:-1],color='black',linewidth=1.5,label="QCD ($H_{T}'$, m$_{T}^{W}$, $p_{T}^{\mathrm{miss}}$  removed)")
plt.hlines(unc0up,bin_edges[:-1],bin_edges_end[:-1],color='black',linewidth=0.7,linestyle='dashed')
plt.hlines(unc0do,bin_edges[:-1],bin_edges_end[:-1],color='black',linewidth=0.7,linestyle='dashed')

#plt.hlines(h1,bin_edges[:-1],bin_edges_end[:-1],color='blue',linewidth=1.5,label=r"QCD ($H_{T}'$, m$_{T}^{W}$  removed)")
plt.hlines(unc1up,bin_edges[:-1],bin_edges_end[:-1],color='blue',linewidth=0.7,linestyle='dashed')
plt.hlines(unc1do,bin_edges[:-1],bin_edges_end[:-1],color='blue',linewidth=0.7,linestyle='dashed')

#plt.hlines(h2,bin_edges[:-1],bin_edges_end[:-1],color='red',linewidth=1.5,label=r"QCD ($H_{T}'$ removed)")
plt.hlines(unc2up,bin_edges[:-1],bin_edges_end[:-1],color='red',linewidth=0.7,linestyle='dashed')
plt.hlines(unc2do,bin_edges[:-1],bin_edges_end[:-1],color='red',linewidth=0.7,linestyle='dashed')


hep.histplot(h0, histtype='step',
    linewidth=1.5, edgecolor='black',alpha=0.7, label=r"QCD ($H_{T}'$, m$_{T}^{W}$, $p_{T}^{\mathrm{miss}}$  removed)")#,density=True)
'''
hep.histplot(unc0up, histtype='step',
    linewidth=1.5, edgecolor='black', linestyle='dashed', alpha=0.7)#,density=True)

hep.histplot(unc0do, histtype='step',
    linewidth=1.5, edgecolor='black', linestyle='dashed', alpha=0.7)#,density=True)
'''
hep.histplot(h1, histtype='step',
    linewidth=1.5, edgecolor='blue',  alpha=0.7, label=r"QCD ($H_{T}'$, m$_{T}^{W}$ removed)")#,density=True)
'''
hep.histplot(unc1up, histtype='step',
    linewidth=1.5, edgecolor='blue',  alpha=0.7,linestyle='dashed')#,density=True)
hep.histplot(unc1do, histtype='step',
    linewidth=1.5, edgecolor='blue',  alpha=0.7,linestyle='dashed')#,density=True)
'''

hep.histplot(h2, histtype='step',
    linewidth=1.5, edgecolor='red',  alpha=0.7, label=r"QCD ($H_{T}'$  removed)")#,density=True)
'''
hep.histplot(unc2up, histtype='step',
    linewidth=1.5, edgecolor='red',  alpha=0.7,linestyle='dashed')#,density=True)
hep.histplot(unc2do, histtype='step',
    linewidth=1.5, edgecolor='red' ,alpha=0.7,linestyle='dashed')#,density=True)
'''

custom_labels = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7']  # Corresponding labels
plt.xticks(bin_edges[:-1], custom_labels)  # Rotate labels if needed



# Customize the plot using mplhep's CMS style
plt.xlabel(f"{var}")  # Label for x-axis, adjust units if applicable
plt.xlabel("MVA score")
plt.ylabel("Events")
#plt.title("Summed Histogram of QCD in CR")

# Add CMS style text
#hep.cms.text("Preliminary", loc=0)
hep.cms.label(llabel='Internal', data=False, lumi=0.302,year=2017,com=5.02)

# Add legend and show the plot
plt.legend()
plt.savefig(f"QCD_forEVAN_AN/{var}_lines.png")
