import pickle
import gzip
import numpy as np
import matplotlib.pyplot as plt

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/cafea/relaxed_cuts_JES10/QCD.pkl.gz','rb') as file: #/mnt_pool/c3_users/user/acerom/github/tchannel-5TeV/MVA_marcos/
	loaded_data=pickle.load(file)
	

channels=['e','m']
levels=['2j0b','2j1b','3j1b','3j2b']
levels=['3j2b']
#levels=['3j1b','4j1b','g5j1b','3j2b','4j2b','g5j2b']

x=np.array([1,2,3,4,5])
for lev in levels:
	
	
	for ch in channels:
		y_up=np.array([])
		y_do=np.array([])
		for var in [2,3,4,5,6]:
			
			if var!=6:upper=loaded_data['counts'].integrate('syst',f'QCDUp_{var}').values()[('QCD',ch,lev)]
			central=loaded_data['counts'].integrate('syst','norm').values()[('QCD',ch,lev)]
			if var!=6:down=loaded_data['counts'].integrate('syst',f'QCDDown_{var}').values()[('QCD',ch,lev)]
			if var==6:upper=loaded_data['counts'].integrate('syst','QCDUp').values()[('QCD',ch,lev)];down=loaded_data['counts'].integrate('syst','QCDDown').values()[('QCD',ch,lev)]
			countup=abs((upper-central)/central*100)
			countdown=abs((down-central)/central*100)
		#	print('upper unc in',ch,lev, 'is:',countup)
		#	print('down unc in',ch,lev, 'is:',countdown)
		#	print('mean unc (value taken to do the plot):',0.5*(countup+countdown))
		#	print(ch,lev,'central',central*302,'upper',upper*302,'down',down*302)
		#	print('\n','\n')
			y_up=np.append(y_up,countup)
			y_do=np.append(y_do,countdown)
			if ch=='e':y_el_up=y_up;y_el_do=y_do;central_e=(central*302).item()
			else:y_m_up=y_up;y_m_do=y_do;central_m=(central*302).item()
		
	fig, ax = plt.subplots()
	ax.scatter(x, y_m_up, color='darkgreen', marker='o', s=25, edgecolors='black', label='Muon Up')
	ax.plot(x, y_m_up, color='darkgreen', linestyle='-', linewidth=1)
	ax.scatter(x, y_m_do, color='limegreen', marker='^', s=25, edgecolors='black', label='Muon Down')
	ax.plot(x, y_m_do, color='limegreen', linestyle='-', linewidth=1)
	ax.scatter(x, y_el_up, color='darkviolet', marker='s', s=25, edgecolors='black', label='Electron Up')
	ax.plot(x, y_el_up, color='darkviolet', linestyle='-', linewidth=1)
	ax.scatter(x, y_el_do, color='orchid', marker='D', s=25, edgecolors='black', label='Electron Down')
	ax.plot(x, y_el_do, color='orchid', linestyle='-', linewidth=1)
	ax.set_title(f'{lev}           Nom. e: {central_e:.1f} Nom. mu: {central_m:.1f}', fontsize=16)
	ax.set_xticks([1,2,3,4,5])
	ax.set_xlabel('N GeV variation', fontsize=14)
	ax.set_ylabel('% QCD variation', fontsize=14)	
	#ax.set_ylim(0,200)	
	
	ax.legend(loc='best', fontsize=12)

# Adding grid lines
	ax.grid(True, which='both', linestyle='--', linewidth=0.5)

# Enhancing the plot aesthetics
	ax.set_facecolor('#f9f9f9')
	fig.patch.set_facecolor('#f9f9f9')	
	plt.savefig('QCD_unc_{}.png'.format(lev))


'''
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
