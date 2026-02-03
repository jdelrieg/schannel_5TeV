import matplotlib.pyplot as plt
import mplhep as hep
import gzip
import pickle

with gzip.open('/mnt_pool/c3_users/user/jriego/tchannel5TeV/splitting_tchan_tbar/prueba/TTPS_prueba.pkl.gz', 'rb') as file:
    loaded_data = pickle.load(file)
    
x=loaded_data['2j1b_jeteta'].value #2j1b_btagsf
y=loaded_data['2j1b_btagscore'].value
z=loaded_data['2j1b_btagsf'].value


plt.style.use(hep.style.CMS)

#plt.scatter(x, y, color='blue', marker='o',s=2)
#plt.ylim(0,0.12)

plt.xlabel(r'$\eta$')
#plt.xlabel(r'$\eta_{\mathrm{jet}}$')
plt.ylabel(r'Score', labelpad=15)
plt.legend()

#plt.hist(y, color='blue', edgecolor='black')
#plt.xlim(0.93,1.05)

plt.hist2d(x, y, bins=([-5,-4,-3,-2.5,-2,0,2,2.5,3,4,5], [-2,0,0.1,0.2,0.3,0.4,0.5,1]), cmap='viridis')

# Add color bar
plt.colorbar(label='Frequency')


plt.savefig("prueba/2dhist.png")
plt.clf()    
