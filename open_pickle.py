import pickle
import gzip
import copy
import numpy as np

def save_histograms(histograms, output_path):
    with gzip.open(output_path, 'wb') as f:
        pickle.dump(histograms, f)

with gzip.open('/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/ht_low/TTPS_part3.pkl.gz','rb') as file: #
	loaded_data=pickle.load(file)

'''
output_file='QCD_2.pkl.gz'
channels=['e_plus','e_minus','m_plus','m_minus']
levels=['2j1b','3j1b','3j2b']
#levels=['3j1b','4j1b','g5j1b','3j2b','4j2b','g5j2b']


# Make a copy of the histograms to avoid modifying the original
modified_hists = copy.deepcopy(loaded_data)

for channel in channels:
	for level in levels:
		modified_hists['MVAscore_relaxed_b10'].values()['QCD',channel,level,'norm'][:10]=np.zeros(10);modified_hists['MVAscore_relaxed_b10'].values()['QCD',channel,level,'norm'][-10:]=modified_hists['MVAscore_relaxed_b10'].values()['QCD',channel,level,'norm'][-10:]

save_histograms(modified_hists, output_file)

for ch in channels:
	for lev in levels:
		upper=loaded_data['counts'].integrate('syst','QCDUp').values()[('QCD',ch,lev)]
		central=loaded_data['counts'].integrate('syst','norm').values()[('QCD',ch,lev)]
		down=loaded_data['counts'].integrate('syst','QCDDown').values()[('QCD',ch,lev)]
		countup=abs((upper-central)/central*100)
		countdown=abs((down-central)/central*100)
		print('upper unc in',ch,lev, 'is:',countup)
		print('down unc in',ch,lev, 'is:',countdown)
		print('mean unc (value taken to do the plot):',0.5*(countup+countdown))
		print(ch,lev,'central',central*302,'upper',upper*302,'down',down*302)
		print('\n','\n')


channels=['e','m']
levels=['2j0b','2j1b','3j1b','3j2b']
#levels=['3j1b','4j1b','g5j1b','3j2b','4j2b','g5j2b']
for ch in channels:
	for lev in levels:
		upper=loaded_data['counts_nheavyjets'].integrate('syst','norm').values()[('W0JetsToLNu',ch,lev)]
		central=loaded_data['counts'].integrate('syst','norm').values()[('W0JetsToLNu',ch,lev)]
		countup=abs(upper/central*100)
		print('percentage of heavy in',ch,lev, 'is:',countup)
		#print('down unc in',ch,lev, 'is:',countdown)
		#print('mean unc (value taken to do the plot):',0.5*(countup+countdown))
		print('\n','\n')
'''



exit()

tt=0
tW=0
tchan=0
dy=0
wjets=0
data=0
suma=0
with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/DY_M10to50.pkl.gz','rb') as file: #DY_M10to50
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('DYJetsToLLM10to50', 'm', '2j1b', 'norm')]
dy=dy+loaded_data['ht'].sum('ht').values()[('DYJetsToLLM10to50', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/DY_M50.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('DYJetsToLLMLL50', 'm', '2j1b', 'norm')]
dy=dy+loaded_data['ht'].sum('ht').values()[('DYJetsToLLMLL50', 'm', '2j1b', 'norm')]


with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tbarchannel_part0.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tbarchannel', 'm', '2j1b', 'norm')]
tchan=tchan+loaded_data['ht'].sum('ht').values()[('tbarchannel', 'm', '2j1b', 'norm')]


with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tbarchannel_part1.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tbarchannel', 'm', '2j1b', 'norm')]
tchan=tchan+loaded_data['ht'].sum('ht').values()[('tbarchannel', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tbarchannel_part2.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tbarchannel', 'm', '2j1b', 'norm')]
tchan=tchan+loaded_data['ht'].sum('ht').values()[('tbarchannel', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tbarchannel_part3.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tbarchannel', 'm', '2j1b', 'norm')]
tchan=tchan+loaded_data['ht'].sum('ht').values()[('tbarchannel', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tbarchannel_part4.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tbarchannel', 'm', '2j1b', 'norm')]
tchan=tchan+loaded_data['ht'].sum('ht').values()[('tbarchannel', 'm', '2j1b', 'norm')]


with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tbarchannel_part5.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tbarchannel', 'm', '2j1b', 'norm')]
tchan=tchan+loaded_data['ht'].sum('ht').values()[('tbarchannel', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tbarchannel_part6.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tbarchannel', 'm', '2j1b', 'norm')]
tchan=tchan+loaded_data['ht'].sum('ht').values()[('tbarchannel', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tbarW.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tbarW', 'm', '2j1b', 'norm')]
tW=tW+loaded_data['ht'].sum('ht').values()[('tbarW', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tchannel_part0.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tchannel', 'm', '2j1b', 'norm')]
tchan=tchan+loaded_data['ht'].sum('ht').values()[('tchannel', 'm', '2j1b', 'norm')]


with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tchannel_part1.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tchannel', 'm', '2j1b', 'norm')]
tchan=tchan+loaded_data['ht'].sum('ht').values()[('tchannel', 'm', '2j1b', 'norm')]


with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tchannel_part2.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tchannel', 'm', '2j1b', 'norm')]
tchan=tchan+loaded_data['ht'].sum('ht').values()[('tchannel', 'm', '2j1b', 'norm')]


with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tchannel_part3.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tchannel', 'm', '2j1b', 'norm')]
tchan=tchan+loaded_data['ht'].sum('ht').values()[('tchannel', 'm', '2j1b', 'norm')]


with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tchannel_part4.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tchannel', 'm', '2j1b', 'norm')]
tchan=tchan+loaded_data['ht'].sum('ht').values()[('tchannel', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tchannel_part5.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tchannel', 'm', '2j1b', 'norm')]
tchan=tchan+loaded_data['ht'].sum('ht').values()[('tchannel', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tchannel_part6.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tchannel', 'm', '2j1b', 'norm')]
tchan=tchan+loaded_data['ht'].sum('ht').values()[('tchannel', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/TTPS_part0.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]
tt=tt+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]


with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/TTPS_part1.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]
tt=tt+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/TTPS_part2.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]
tt=tt+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/TTPS_part3.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]
tt=tt+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/TTPS_part4.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]
tt=tt+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/TTPS_part5.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]
tt=tt+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/TTPS_part6.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]
tt=tt+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/TTPS_part7.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]
tt=tt+loaded_data['ht'].sum('ht').values()[('ttPS', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tW_part0.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tW', 'm', '2j1b', 'norm')]
tW=tW+loaded_data['ht'].sum('ht').values()[('tW', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tW_part1.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tW', 'm', '2j1b', 'norm')]
tW=tW+loaded_data['ht'].sum('ht').values()[('tW', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tW_part2.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tW', 'm', '2j1b', 'norm')]
tW=tW+loaded_data['ht'].sum('ht').values()[('tW', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/tW_part3.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('tW', 'm', '2j1b', 'norm')]
tW=tW+loaded_data['ht'].sum('ht').values()[('tW', 'm', '2j1b', 'norm')]


with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/W0JetsToLNu.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('W0JetsToLNu', 'm', '2j1b', 'norm')]
wjets=wjets+loaded_data['ht'].sum('ht').values()[('W0JetsToLNu', 'm', '2j1b', 'norm')]


with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/W1JetsToLNu.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('W1JetsToLNu', 'm', '2j1b', 'norm')]
wjets=wjets+loaded_data['ht'].sum('ht').values()[('W1JetsToLNu', 'm', '2j1b', 'norm')]


with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/W2JetsToLNu.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('W2JetsToLNu', 'm', '2j1b', 'norm')]
wjets=wjets+loaded_data['ht'].sum('ht').values()[('W2JetsToLNu', 'm', '2j1b', 'norm')]


with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/W3JetsToLNu.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('W3JetsToLNu', 'm', '2j1b', 'norm')]
wjets=wjets+loaded_data['ht'].sum('ht').values()[('W3JetsToLNu', 'm', '2j1b', 'norm')]


with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/QCD.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

suma=suma+loaded_data['ht'].sum('ht').values()[('QCD', 'm', '2j1b', 'norm')]
qcd=(loaded_data['ht'].sum('ht').values()[('QCD', 'm', '2j1b', 'norm')])

#data
with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/HighEGJet.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

data=data+loaded_data['ht'].sum('ht').values()[('HighEGJet', 'm', '2j1b', 'norm')]

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/nocuts_wjets_split/SingleMuon.pkl.gz','rb') as file:
	loaded_data=pickle.load(file)

data=data+loaded_data['ht'].sum('ht').values()[('SingleMuon', 'm', '2j1b', 'norm')]


print('QCD',qcd*302)
print('total',suma*302)
print('tt',tt*302)
print('tW',tW*302)
print('tchan',tchan*302)
print('DY',dy*302)
print('Wjets',wjets*302)
print('data',data)
