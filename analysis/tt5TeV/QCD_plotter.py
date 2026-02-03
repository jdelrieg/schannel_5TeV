import pickle
import gzip
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
from coffea import hist, processor

with gzip.open('/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/Questions/QCD/QCD_noht_nomtw_nomet_cuts.pkl.gz','rb') as file: #/mnt_pool/c3_users/user/acerom/github/tchannel-5TeV/MVA_marcos/
	loaded_data=pickle.load(file)

with gzip.open('/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/split_charge_goodJECs_mistag_comb/QCD_shape/QCD_shapes.pkl.gz','rb') as file: #
	loaded_data=pickle.load(file)
	
hep.style.use("CMS")

var='j0pt'
channels=['e_fake_minus','e_fake_plus','m_fake_minus','m_fake_plus']
levels=['2j1b','3j1b','3j2b']

h0=np.zeros(len(loaded_data[var].values()[('QCD', 'm_fake_minus', '2j1b', 'norm')]))


for ch in channels:
	for level in levels:
		if ch in ['e_fake_minus','e_fake_plus']:
			if (ch=='e_fake_plus')&(level=='3j2b'):
				print('canla') 
			else:
				 h0+=loaded_data[var].values()[('QCD', ch, level, 'norm')]
		else:
			h0+=loaded_data[var].values()[('QCD', ch, level, 'norm')] 

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


# Plot the histogram of h0 with some customization
hep.histplot(h0, histtype='fill',
    linewidth=1.5, edgecolor='black', facecolor='skyblue', alpha=0.7, label="QCD (summed)")

# Customize the plot using mplhep's CMS style
plt.xlabel(f"{var} [GeV]")  # Label for x-axis, adjust units if applicable
plt.ylabel("Events")
#plt.title("Summed Histogram of QCD in CR")

# Add CMS style text
hep.cms.text("Preliminary", loc=0)
hep.cms.lumitext("5.02 TeV")

# Add legend and show the plot
plt.legend()
plt.savefig(f"Questions/QCD/{var}_tricked.png")
