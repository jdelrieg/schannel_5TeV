import pickle
import gzip
import copy
from coffea.hist import Hist
import numpy as np

# Function to load the .pkl.gz file
def load_histograms(file_path):
    with gzip.open(file_path, 'rb') as f:
        return pickle.load(f)

# Function to save the modified histograms to a new .pkl.gz file
def save_histograms(histograms, output_path):
    with gzip.open(output_path, 'wb') as f:
        pickle.dump(histograms, f)


########################################## Initialization of factors: driven from Save_QCD. To be incorporated automatized
fake_rates={('e_minus','2j1b'):0.901,
	('e_plus','2j1b'):0.901,
	('m_minus','2j1b'):0.048,
	('m_plus','2j1b'):0.048,
	('e_minus','3j1b'):0.901,
	('e_plus','3j1b'):0.901,
	('m_minus','3j1b'):0.048,
	('m_plus','3j1b'):0.048,
	('e_minus','3j2b'):0.901,
	('e_plus','3j2b'):0.901,
	('m_minus','3j2b'):0.048,								#These are the Fake rates (parte denominador)
	('m_plus','3j2b'):0.048,
		}

cr_rates={('e_minus','2j1b'):25.317,
	('e_plus','2j1b'):21.648,
	('m_minus','2j1b'):154.044,
	('m_plus','2j1b'):164.443,
	('e_minus','3j1b'):24.756,
	('e_plus','3j1b'):45.339,
	('m_minus','3j1b'):141.927,
	('m_plus','3j1b'):166.839,								#These are the normalizations in each category (N_dat^{CR}-N_MC^CR). In yields or counts
	('e_minus','3j2b'):6.1298,
	('e_plus','3j2b'):3.827,
	('m_minus','3j2b'):11.4098,
	('m_plus','3j2b'):8.183,
		}

joint_rates={key: fake_rates[key] * cr_rates[key]/302 for key in fake_rates}   #operate both factors  AND divide by lumi. This could be done (to be more rigid in the ordering of operations) when the 'CR factors' are calculated

##############################################################################			Up to here the factors

input_file="split_charge_goodJECs_mistag_comb_0QCDlepton/QCD_shape/QCD_shapes.pkl.gz"            #point to input file (driven from QCD mc after passing the same process as other samples, only difference the relaxation of cuts)
output_file = "split_charge_goodJECs_mistag_comb_0QCDlepton/binnedQCD/QCD_ll.pkl.gz"

# Load the histograms
hists = load_histograms(input_file)

# Make a copy of the histograms to avoid modifying the original
modified_hists = copy.deepcopy(hists)

##########################################################################				First part: integrate all channels and levels to get a single shape for all the categories. 
ch=['e_fake_minus','e_fake_plus','m_fake_minus','m_fake_plus']						 # we are integrating in fake channels since those are the enriched by the selection
ch_el=['e_fake_minus','e_fake_plus']
ch_mu=['m_fake_minus','m_fake_plus']

lv=['2j1b','3j1b','3j2b']

for name, hist in modified_hists.items():
	if isinstance(hist, Hist):  # Ensure it's a coffea Hist object
		if not hist.values():
			continue          #Skip if we find not filled histograms (like dummy)
			
		
		if name in ['MVAtth','MVAwp','pttrig','etatrig','abseeta','absmeta']:#'ept','eeta','mpt','meta',
			continue         # for the moment we skip the 'e' only or 'mu' only (lep 0 will serve for this)
			
		if name not in ['ept','eeta']: h0=np.zeros(len(hist.values()[('QCD', 'm_fake_minus', '2j1b', 'norm')])) #initialize to 0 the array
		else:h0=np.zeros(len(hist.values()[('QCD', 'e_fake_minus', '2j1b', 'norm')])) 
		
		for c in ch:
			for l in lv:
				key=('QCD', c, l, 'norm')          #loop through channels and levels (and define the key meanwhile) to sum the entries in each category
				if key not in hist.values():
					continue                       #check if the key wasn't filled while constructing the QCD_original, skip that one and keep summing the filled ones. If this is not added, an error will raise
				h0+=hist.values()[key]
		#once this loop is finished, we have the UNIQUE (integrated over the 12 cats) SHAPE that has to be set as starting point in all the 12 categories

############################################# exception for e only and mu only		
		
		if name in ['mpt','meta']:
			for c1 in ch_el:
				for l in lv:
					key=('QCD', c1, l, 'norm')
					#if key not in hist.values():
					#	continue
					if name=='mpt':h0+=modified_hists['ept'].values()[key];
					else: h0+=modified_hists['eeta'].values()[key]
					
		if name in ['ept','eeta']:
			for c1 in ch_mu:
				for l in lv:
					key=('QCD', c1, l, 'norm')
					#if key not in hist.values():
					#	continue
					if name=='ept':h0+=modified_hists['mpt'].values()[key]
					else: h0+=modified_hists['meta'].values()[key]
################################################################################ end of exception for e only and mu only	
		
		#In the coming loop, we will go through the signal categories and assign to all of them the UNIQUE SHAPE just calculated
		for (channel, level), fake_rate in joint_rates.items(): 
			key=('QCD', channel, level, 'norm')
			
			if key not in hist.values():	                        #From here
				fill_args = {}
				for axis in hist.axes():
					if axis.name == "sample":
						fill_args["sample"] = "QCD"
					elif axis.name == "channel":
						fill_args["channel"] = channel
					elif axis.name == "level":
						fill_args["level"] = level
					elif axis.name == "syst":
						fill_args["syst"] = "norm"
					elif axis.name == name: 
						fill_args[axis.name] = np.zeros(1)
				hist.fill(**fill_args)                              #To here we are just tricking if in the original signal categories (e.g. e_plus 2j1b) the histogram wasn't filled
																    # in that case, we create the histogram for that category and fill it with dummy values
				
			hist.values()['QCD',channel,level,'norm'][:]=h0       #here is where we assign the Unique Shape to all the categories (both the ones which were already filled with very low values and the empty ones that now contain dummy values)
																  #VERY IMPORTANT is the [:]. If not in coffea the values of a histogram dont get updated
		

################################################################################## end of first part. At this point we have a Hist with a unique shape per varaible in all the categories


################################################################################## Second part: the scalings by factors 

for name, hist in modified_hists.items():
	if isinstance(hist, Hist):  # Ensure it's a coffea Hist object
		
		if 'channel' in hist.axes() and 'level' in hist.axes():
			for (channel, level), fake_rate in joint_rates.items():
				
				valores=hist.integrate('channel',channel).integrate('level',level).values() #Do the integration of the shape...
				if valores !={}:
					integral=np.sum(valores[('QCD','norm')]) #... if we have smth, if not...
				else:
					integral=1     #   ... a naive 1
					
				hist.scale({(channel,level):fake_rate/integral},axis=('channel','level'))  #The actual scaling by the rates and division by integral. NOTE that the integral has to be calculated before the scaling
				
				

# Save the modified histograms to a new file
save_histograms(modified_hists, output_file)

print(f"Modified histograms saved to {output_file}")
