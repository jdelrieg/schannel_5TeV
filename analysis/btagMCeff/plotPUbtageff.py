import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import pickle
import gzip

path='/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/histos/btagMCeff.pkl.gz'
outpath='prueba/btag/94X/'
with gzip.open(path,'rb') as file: #
	loaded_data=pickle.load(file)
'''
def DrawPUEff(path, WP='medium', flav='b', year=None, outpath=''):
  var='pileup'
  num=loaded_data[var].values()[(WP, flav)]
  den=loaded_data[var].values()[('all', flav)]
  ratio=num/den
  unc=unc=(np.sqrt(num)*den+num*np.sqrt(den))/(den**2)

  edges=loaded_data[var].axes()[2].edges()
  midpoints=(edges[:-1]+edges[1:])/2
  
  plt.style.use(hep.style.CMS)
  plt.errorbar(midpoints, ratio, yerr=unc, fmt='o', label='Data with uncertainties', color='blue', ecolor='red', capsize=5)
  
  plt.xlabel(r'$N_{vtx}$')
  plt.ylabel('B-tag efficiency') 
  plt.grid(True)
  plt.legend()
  hep.cms.label(llabel='Internal', data=False,year=2017,com=13)
  plt.savefig(f"{outpath}BtagEff_{var}{WP}{flav}.png")

 

for wp in ['loose', 'medium', 'tight']:
  for f in ['b', 'l']:
    DrawPUEff(path, wp, f, outpath=outpath)
'''


def DrawPUEff(path, WP='medium', flav='b', year=None, outpath=''):
    var = 'pileup'
    
    # Get numerator and denominator
    num = loaded_data[var].values()[(WP, flav)]
    den = loaded_data[var].values()[('all', flav)]
    
    # Calculate ratio and uncertainties
    ratio = num / den
    unc = (np.sqrt(num) * den + num * np.sqrt(den)) / (den**2)

    # Get bin edges and calculate midpoints
    edges = loaded_data[var].axes()[2].edges()
    midpoints = (edges[:-1] + edges[1:]) / 2

    return midpoints, ratio, unc

def plot_btag_eff(path, flav, outpath=''):
    # Define colors for each WP
    wp_colors = {
        'loose': 'blue',
        'medium': 'green',
        'tight': 'red'
    }
    
    plt.figure(figsize=(6, 6))
    
    for wp in ['loose', 'medium', 'tight']:
        #if wp != 'medium': continue
        midpoints, ratio, unc = DrawPUEff(path, WP=wp, flav=flav)
        plt.errorbar(midpoints, ratio, yerr=unc, fmt='o', label=f'WP={wp}', color=wp_colors[wp], ecolor=wp_colors[wp], capsize=5)
    
    plt.xlabel(r'$N_{vtx}$')
    plt.ylabel('B-tag efficiency') if flav=='b' else plt.ylabel('udsg misId')
    if flav=='c': plt.ylabel('c misId');plt.text(0.05, 0.5, r'Mis Id rate in $\sqrt{s}=$5.02 TeV: 0.094', transform=plt.gca().transAxes, fontsize=12,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')) 
    if flav=='l': plt.yscale('log')
    #plt.ylim(-0.02,0.35) if flav=='l' else plt.ylim(0.2,1.1)
    plt.grid(True)
    plt.legend()
    
    # Add CMS label
    hep.cms.label(llabel='Internal', data=False, year=2017, com=13)

    #plt.text(0.05, 0.85, r'B tag Eff. in $\sqrt{s}=$5.02 TeV: 0.741', transform=plt.gca().transAxes, fontsize=12,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')) if flav=='b' else plt.text(0.08, 0.75, r'Mis Id rate in $\sqrt{s}=$5.02 TeV: 0.0047', transform=plt.gca().transAxes, fontsize=12,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

    
    # Save as a separate file for each flavor
    plt.tight_layout()
    plt.savefig(f"{outpath}BtagEff_{flav}.png")
    plt.close()

# Generate and save separate plots for b and l flavors
plot_btag_eff(path, 'b', outpath=outpath)
plot_btag_eff(path, 'l', outpath=outpath)
plot_btag_eff(path, 'c', outpath=outpath)