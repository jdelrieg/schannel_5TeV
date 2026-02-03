'''
 Same as ControlPlots.py, but avoids loading QCD for drawing non-QCD iso and non-iso plots
'''


from config import *
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
from cafea.plotter.plotter import DrawUncPerBin

from ControlPlots import DrawAll

def QCDnormUncPlot(path, outpath='./'):
   qcd = GetHisto(path, 'counts', {'sample':'QCD'})
   levels = ['2j1b']#['3j1b', '3j2b', '4j1b', '4j2b', 'g5j1b', 'g5j2b']
   channels = ['e', 'm']
   syst = 'QCD'
   for c in channels:
    for l in levels:
      hnorm = qcd.integrate('channel', c).integrate('level', l).integrate('syst', 'norm')
      if hnorm.values() == {}:
        maxrel = 0
      else:
        norm = hnorm.values()[()]
        up   = qcd.integrate('channel', c).integrate('level', l).integrate('syst', syst+'Up').values()[()]
        do   = qcd.integrate('channel', c).integrate('level', l).integrate('syst', syst+'Down').values()[()]
        maxrel = max(abs(up-norm), abs(do-norm))/norm
      print('QCD norm unc in %s %s: %.2f'%(c, l, maxrel))

if __name__=="__main__":
  outpathQCD = 'selection_mtw_reliso/Control_mas_QCD/'
  if not os.path.isdir(outpathQCD): os.makedirs(outpathQCD)

  #QCDnormUncPlot(path+'/QCD.pkl.gz', outpathQCD)
  #exit()

  print(' >> Loading histograms! This might take a while...')
  plt = plotter(path, prDic=processDic_noQCD, bkgList=bkglist_noQCD, colors=colordic, lumi=lumi)
  plt.SetLumi(lumi, "pb$^{-1}$", "5.02 TeV")
  plt.SetRatio(True)
  plt.SetVerbose(verbose)
  plt.plotData = True

  variables = ['mtwnocut','metnocut']#['mtw','met', 'ht', 'mt', 'j0eta', 'j0pt', 'ept', 'eeta', 'mpt', 'meta','ht_atlas','mlb','met_mtw','metnocut']
  levels = ['2j1b']#['3j1b', '3j2b', '4j1b', '4j2b', 'g5j1b', 'g5j2b']
  channels_fake = ['e_fake', 'm_fake']
  channels = ['e', 'm']
  #channels_fake=['lep_pluss_fake','lep_minus_fake']
  #channels=['lep_pluss','lep_minus']



  ### Uncertainties
  #outpath = outpathQCD+'/uncertainty/'
  #if not os.path.isdir(outpath): os.makedirs(outpath)
  #outname = 'QCD_uncertainty'
  #fname = path + 'QCD.pkl.gz'
  #DrawUncPerBin(fname,  ['QCD'], process='QCD', var='ht', outpath=outpath, outname=outname, prDic=processDic, cat={'channel':'e', 'level':'3j1b'})
  #exit()

  ### Non-iso plots
  outpath = outpathQCD+'/non-iso/'
  if not os.path.isdir(outpath): os.makedirs(outpath)
  print('[No-iso control plots] Output = ', outpath)
  DrawAll(plt, outpath, variables, levels, channels_fake, nSlots=nSlots)

  ### No QCD plots
  outpath = outpathQCD+'/noQCD/'
  if not os.path.isdir(outpath): os.makedirs(outpath)
  print('[No-QCD control plots] Output = ', outpath)
  DrawAll(plt, outpath, variables, levels, channels, nSlots=nSlots)


