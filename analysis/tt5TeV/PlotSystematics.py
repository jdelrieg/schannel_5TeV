'''
 Usage: python ljets/plotSyst.py -p histos/2jun2022_btagM_jet25/TT.pkl.gz
'''

from config import *
#datatoday='24mar23'

#baseweb='/nfs/fanae/user/jriego/www/public/tt5TeV/'
outpath = baseweb+datatoday+'/systematics/'
outpath='systematics_tchannel/'

print(' >> Output = ', outpath)
if not os.path.isdir(outpath): os.makedirs(outpath)

print(' >> Loading histograms! This might take a while...')
plt = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic,outpath=outpath, lumi=lumi)
plt.SetLumi(lumi, "pb$^{-1}$", "5.02 TeV")
plt.SetRatio(True)
plt.SetYRatioTit('Ratio')
plt.plotData = True
plt.SetOutput(output)

variables = ['counts','met', 'medianDRjj', 'MVAscore', 'ht', 'st', 'counts', 'njets', 'nbtags', 'met', 'j0pt', 'j0eta', 'ept', 'eeta', 'mpt', 'meta','mjj', 'mt', 'ptjj', 'minDRjj', 'medianDRjj', 'u0pt', 'u0eta', 'minDRuu', 'medianDRuu', 'ptlb', 'ptuu', 'mlb', 'sumallpt', 'dRlb']
#levels = ['g4jets', '3j1b', '3j2b', '4j1b', '4j2b', 'g5j1b', 'g5j2b']
levels = ['2j1b','3j1b','3j2b']

channels = ['e', 'm']
#systematics =  ['FSR', 'ISR', 'JES', 'btagSF', 'lepSF', 'prefire']
systematics = ['ISR','FSR', 'trigSF', 'btagSF', 'elecSF', 'muonSF', 'prefire','JER','MC', 'AbsStat', 'AbsScale', 'AbsMPF', 'Frag', 'ECAL', 'HCAL', 'Flavor', 'RelStat', 'RelPt', 'RelBal', 'RelJER', 'L3Res','MET_UnclusteredEnergy']

plt.SetSystematics(systematics)


def DrawAllForVar(var, level=None):
  xtit = RebinVar(plt, var, level)
  for l in levels:
    for c in channels:
      for s in systematics:
        outname = '%s_%s_%s_%s_%s'%(pr, var, c, l, s)
        plt.SetOutput(outname)
        DrawSyst(var, s, pr, c, l)
 

def DrawSyst(var, syst, process='tt', chan='m', level='g4jets'): 
  colors = ['k', 'r', 'b']
  labels = ['Nominal', '%s up'%syst, '%s down'%syst]
  h = plt.GetHistogram(var)
  print('h = ', h)
  PrintHisto(h)
  systlist = [x.name for x in h.identifiers('syst')]
  upsyst = syst+'Up'; dosyst = syst+'Down' if syst+'Down' in systlist else syst+'Do';
  if not upsyst in systlist or not dosyst in systlist:
    print (" >> WARNING: No syst %s found -- List: "%(syst), systlist)
    return
  selec = [{'channel':chan, 'level':level, 'syst':'norm'}, {'channel':chan, 'level':level, 'syst':upsyst}, {'channel':chan, 'level':level, 'syst':dosyst} ]
  plt.SetRatioRange(ymin=0.90, ymax=1.10)
  plt.DrawComparison(var, process, selec, labels, colors, [], [])

  


pr = 'tt'
pr='tchan'
if not var is None:
  if systch == 'hdamp':
    hdampup,hdampdo = GetModSystHistos(path, 'TT_hdamp', 'hdamp', var=var)
    plt.AddExtraBkgHist([hdampup, hdampdo], add=True)
  elif systch == 'UE':
    tuneup , tunedo = GetModSystHistos(path, 'TT_UE', 'UE', var=var)
    plt.AddExtraBkgHist([tuneup, tunedo], add=True)
  DrawSyst(var, systch, pr, ch, level)
  exit()

else:
  var = ['MVAscore']#['met']
  for v in var:
    DrawAllForVar(v)
  exit()

