from config import *
from cafea.plotter.plotter import plotter

### Values to be updated w.r.t. TOP-20-004
tWunc = 0.2 # This one is correlated with l+jets
lumiUnc = 0.015

#### Values from TOP-20-004
bkg = ['tW', 'Nonprompt', 'DY', 'VV']
signal = 'tt'
yields = {'tt':187., 'tW':8., 'Nonprompt':2., 'DY':10., 'VV':4., 'data':194}

bkgunc = {'tW':tWunc, 'Nonprompt':0.3, 'DY':0.3, 'VV':0.3}

expUnc = {'muonSF':0.006, 'eleSF':0.016, 'trigSF':0.013, 'prefire':0.014, 'JES':0.022, 'JER':0.012}
modelingUnc = {'FSR':0.011, 'ISR':0.001, 'PDF':0.003, 'Scales':0.002, 'hdamp':0.010, 'UE':0.007} 

outpath = path + '/fromTOP20004/'
if not os.path.exists(outpath):
  os.makedirs(outpath)
outname = 'counts_em_g2jets'

def CreateCoffeaHist(outname='counts_em_g2jets'):
  ''' Save TOP-20-004 yields and "shape" uncertainties in a coffea histogram '''
  from coffea import hist
  
  ### Create the histogram
  h = hist.Hist("Events", hist.Cat("process", "process"), hist.Cat('syst', 'syst'), hist.Cat('channel', 'channel'), hist.Cat('level', 'level'), hist.Bin("counts", "counts", 1, 0., 2.))
  
  ### Fill signal yields
  y_signal = np.array([yields[signal]/1e6]*int(1e6))
  h.fill(**{'syst':'norm', 'weight':y_signal, 'process':signal, 'channel':'em', 'level':'g2jets', 'counts':np.ones_like(y_signal)})
  
  ### Fill data yields
  y_data = np.ones(yields['data'], dtype=float)
  h.fill(**{'syst':'norm', 'weight':y_data, 'process':'data', 'channel':'em', 'level':'g2jets', 'counts':np.ones_like(y_data)})
  
  ### Fill background yields
  for pr in bkg:
    y = np.array([yields[pr]])
    h.fill(**{'syst':'norm', 'weight':y, 'process':pr, 'channel':'em', 'level':'g2jets', 'counts':np.ones_like(y)})
  
  ### Fill signal shape uncertainties (experimental)
  for unc in expUnc:
    y = yields[signal]
    y_up = y*(1+expUnc[unc])
    y_do = y*(1-expUnc[unc])
    y = np.array([yields[signal]*expUnc[unc]])
    h.fill(**{'syst':unc+'Up',   'weight':np.array([y_up]), 'process':signal, 'channel':'em', 'level':'g2jets', 'counts':np.ones_like(np.array([y_up]))})
    h.fill(**{'syst':unc+'Down', 'weight':np.array([y_do]), 'process':signal, 'channel':'em', 'level':'g2jets', 'counts':np.ones_like(np.array([y_do]))})
  
  
  saveHistos(outpath, outname, {'counts':h}, verbose=True)


if __name__ == '__main__':
  # Create in coffea format
  ####################################################################
  fname = outpath + outname + '.pkl.gz'
  if not os.path.exists(fname) or force:
    if verbose: print('Creating coffea file: ', fname)
    CreateCoffeaHist(outname)

  # Create root file for combine
  ####################################################################
  pathcomb = "/nfs/fanae/user/juanr/CMSSW_10_2_13/src/combinescripts/tt5TeVljets/"+datatoday + '/'
  if not os.path.exists(pathcomb):
    os.makedirs(pathcomb)
  outname = "counts_em_g2jets"
  fcombine = pathcomb + outname
  if not os.path.exists(f"{pathcomb+outname}") or force:
    plt = plotter(fname, prDic={}, lumi=1., var='counts')
    plt.SetOutpath(pathcomb)
    plt.SetVerbose(verbose)
    plt.SaveCombine('counts', outname, categories={'level':'g2jets', 'channel':'em'})
    if verbose: print('Created root file for combine: ', fcombine+'.root')

  # Create datacard
  ####################################################################
  from cafea.modules.CreateDatacardFromRootfile import Datacard
  from cafea.plotter.plotter import GetHisto

  # Paths and file names...
  fcombine = fcombine+'.root'
  if  not '/' in fcombine: outpath = './'
  else: outpath = fcombine[:fcombine.rfind('/')]
  oname = fcombine[fcombine.rfind('/')+1:] if '/' in fcombine else fcombine
  if oname.endswith('.root'): oname = oname[:-5]
  if '/' in oname: oname[oname.rfind('/')+1:]
  oname = 'datacard_'+oname
  if not oname.endswith('.txt'): oname += '.txt'


  norm = [bkgunc[pr] for pr in bkg]
  systList = expUnc.keys()
  d = Datacard(fcombine, signal, bkg, lumiUnc, norm, systList, nSpaces=12, verbose=verbose)

  # Add modeling uncertainties
  for mod in modelingUnc:
    d.AddExtraUnc(mod, modelingUnc[mod], signal)  
  d.SetOutPath(outpath)
  if verbose: print('Saving datacard to: %s'%outpath)
  d.Save(oname)
