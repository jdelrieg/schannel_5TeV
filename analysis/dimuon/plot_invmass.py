from cafea.plotter.plotter import *
import matplotlib.pyplot as plt
path = 'histos/DoubleMuon.pkl.gz'

var = 'invmassg'
h = GetHisto(path, var)

# Normalize each bin to 1 GeV
bins, vals = GetXYfromH1D(h, mode='edges', errors=False)
vals = vals / (bins[1:] - bins[:-1])
h = GetH1DfromXY(bins, vals)

fig, ax = plt.subplots(1, 1, figsize=(12,8))
hist.plot1d(h, ax=ax, fill_opts=None, error_opts=None, line_opts={'color':'black', 'linewidth':2}, legend_opts={'loc':'upper right', 'fontsize':20})
ax.set_xlabel('m$_{\mu\mu}$ (GeV)', fontsize=20)
ax.set_ylabel('Events / GeV', fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax.tick_params(axis='both', which='minor', labelsize=20)
#ax.set_title('Dimuon invariant mass spectrum at 5.02 TeV in pp collisions', fontsize=20)

# do not show the legend
ax.legend().set_visible(False)

# Add labels
ax.text(0.02, 0.99, 'CMS', transform=ax.transAxes, fontsize=20, fontweight='bold', va='top')
# Lumi label
ax.text(0.68, 0.99, '296.1 pb$^{-1}$ (5.02 TeV)', transform=ax.transAxes, fontsize=20, va='top')

# labels for the resonances
ax.text(0.009, 0.701, '$\eta$', transform=ax.transAxes, fontsize=20, va='top')
ax.text(0.058, 0.74, r'$\rho,\omega$', transform=ax.transAxes, fontsize=20, va='top')
ax.text(0.122, 0.771, '$\phi$', transform=ax.transAxes, fontsize=20, va='top')
ax.text(0.314, 0.87, '$J/\psi$', transform=ax.transAxes, fontsize=20, va='top')
ax.text(0.355, 0.711, '$\Psi$', transform=ax.transAxes, fontsize=20, va='top')
ax.text(0.527, 0.76, r'$\Upsilon$', transform=ax.transAxes, fontsize=20, va='top')
ax.text(0.941, 0.51, '$Z$', transform=ax.transAxes, fontsize=20, va='top')

if var == 'invmass' or var == 'invmassg':
  ax.set_yscale("log")
  ax.set_xscale("log")
  ax.set_ylim(1, 1e9)
  ax.set_xlim(0.5, 120)

elif var == 'Z':
  ax.set_xlim(80, 100)

plt.tight_layout()
fig.savefig(var+".png")
