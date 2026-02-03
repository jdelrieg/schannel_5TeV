from config import *
from QCD import *
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

plt = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi, var=var)
plt.SetLumi(lumi, "pb$^{-1}$", "5.02 TeV")

def Draw(var, categories, output=None, label='', outpath='temp/', doQCD=False, doRatio=True):
  doQCD = True
  if not CheckHistoCategories(plt.hists[var], categories):
    print(f'Nop... [{var}] cat = ', categories)
    return
  plt.ResetExtraBkg()
  plt.SetRatio(doRatio)
  plt.plotData = doData
  plt.SetCategories(categories)
  label = GetChLab(categories['channel']) + GetLevLab(categories['level']) 
  plt.SetRegion(label)
  plt.SetOutpath(outpath)
  plt.SetOutput(output)
  #plt.SetLogY()
  rebin = None; xtit = ''; ytit = ''
  if doQCD: 
    qcd = QCD(pathQCD, prDic=processDic, bkglist=bkglist, lumi=lumi, categories=categories, var=var)
    hqcd = qcd.GetQCD(var, categories, 0)
    #hqcdup = qcd.GetQCD(var, categories,  1)
    #hqcddo = qcd.GetQCD(var, categories, -1)

  b0 = None; bN = None; binRebin = None
  if var in ['minDRjj', 'minDRuu']:
    b0 = 0.4; bN = 2.0
  elif var in ['medianDRjj']:
    b0 = 1; bN = 4.0
  elif var == "njets" and categories['level'] != 'incl':
    b0 = 4; bN = 10
  #elif var in ['ht']:
  #  b0 = 2
  elif var in ['st']:
    b0 = 120; bN = 600;
  elif var in ['sumallpt']:
    b0 = 0; bN = 200
    xtit = '$\sum_\mathrm{j,\ell}\,\mathrm{p}_{T}$ (GeV)'
  elif var in ['met','u0pt', 'ptuu', 'ptjj', 'metnocut']:
    b0 = 2;
  elif var in ['MVAscore']:
    b0 = 0.1; bN = 0.8
    binRebin = 2
  elif var in ['ht']:
    b0 = 4;

  if b0 is not None:
    plt.SetRebin(var, b0, bN, includeLower=True, includeUpper=True, binRebin=binRebin)
    if doQCD: 
      hqcd   = Rebin(hqcd, var, b0, bN, includeLower=True, includeUpper=True, binRebin=binRebin)
      #hqcdup = Rebin(hqcdup, var, b0, bN, includeLower=True, includeUpper=True)
      #hqcddo = Rebin(hqcddo, var, b0, bN, includeLower=True, includeUpper=True)
    
  if doQCD: 
    #hqcd.add(hqcdup)
    #hqcd.add(hqcddo)
    #print('QCDUp   = ', hqcdup.integrate('process', 'QCD').integrate('syst', 'QCDUp'  ).values(overflow='all'))
    #print('QCDDown = ', hqcddo.integrate('process', 'QCD').integrate('syst', 'QCDDown').values(overflow='all'))
    #print('QCD     = ', hqcd  .integrate('process', 'QCD').integrate('syst', 'norm'   ).values(overflow='all'))
    plt.AddExtraBkgHist(hqcd)

  aname = None
  if   var == 'l0pt': aname = 'lep0pt'
  elif var == 'l0eta': aname = 'lep0eta'
  plt.Stack(var, xtit=xtit, ytit=ytit, aname=aname, dosyst=True)
  #plt.PrintYields('counts')

def Print2lplots():
  for c in ['em', 'ee', 'mm']:
    for l in ['incl', 'g2jets']:
      outp = outpath+'/2l/'+c+'/'+l+'/'
      cat = {'channel':c, 'level':l, 'syst':'norm'}
      for var in ['l0pt','counts', 'njets', 'nbtags', 'ht', 'met', 'j0pt', 'j0eta', 'l0eta', 'invmass', 'invmass2']:
        if l=='incl' and var in ['j0pt', 'j0eta']: continue
        outname = "%s_%s_%s"%(var, c, l)
        Draw(var, cat, outname, outpath=outp)

def Print1lplots(channels, levels):
  if not isinstance(channels, list): channels = [channels]
  if not isinstance(levels, list): levels = [levels]
  outp = outpath+'/1l/'
  for c in channels: #['m', 'e']:#, ['e','m']]: #, 'e_fake', 'm_fake']:
    doQCD = not 'fake' in c
    doRatio = not 'fake' in c
    #for l in ['incl', 'g2jets', 'g4jets', '0b', '1b', '2b', '2j1b', '3j1b', '3j2b', '4j1b', '4j2b', 'g5j1b', 'g5j2b']:
    for l in levels: #['2j1b']:#, '4j1b', '4j2b', 'g5j1b', 'g5j2b']:
      cat = {'channel':c, 'level':l}#, 'syst':'norm'}
      clab = c if not isinstance(c, list) else 'l'
      outp = outpath+'/1l/'+clab+'/'+l+'/'
      for var in ['met']:#['medianDRjj', 'DNNscore', 'ht', 'st', 'counts', 'njets', 'nbtags', 'met', 'j0pt', 'j0eta', 'ept', 'eeta', 'mpt', 'meta','mjj', 'mt', 'ptjj', 'minDRjj', 'medianDRjj', 'u0pt', 'u0eta', 'minDRuu', 'medianDRuu', 'ptlb', 'ptuu', 'mlb', 'sumallpt', 'dRlb']:#, 'DNNscore']: dRlb
        if l=='incl' and var in ['j0pt', 'j0eta']: continue
        if var == 'DNNscore' and not l in ['2j1b', '3j1b', '3j2b']: continue
        outname = "%s_%s_%s"%(var, clab, l)
        Draw(var, cat, outname, outpath=outp, doQCD=doQCD)


if not var is None:
  categories = { 'channel': ch, 'level' : level}#, 'syst':'norm'}
  Draw(var, categories, output, doQCD=True if ((len(ch) <= 1 and str(*ch) in ['e', 'm']) or (ch[0] in ['e', 'm']) ) else False)


else:
  outpatho = '28nov2022/'
  outpath = '/nfs/fanae/user/juanr/www/public/tt5TeV/ljets/' + outpatho
  if not os.path.exists(outpath): os.makedirs(outpath)
  #Print2lplots()
  Print1lplots('m', ['2j1b'])


