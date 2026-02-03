### Usage:
# python analysis/ttLepEff/stack.py -p histos/ttLepEff/ --data
# python analysis/ttLepEff/stack.py -p histos/ttLepEff/ -l incl -c em -v njets --data


from config import *
from DrawEff import GetNonprompt

def Draw(var, categories, output=None, label='', outpath='temp/', doRatio=True, nonpromptDD=False):
  bkglist = ['tt', 'tW', 'DY', 'Diboson', 'WJets', 'semilep']
  if nonpromptDD:
    bkglist = ['tt', 'tW', 'DY', 'Diboson']
  plt = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi, var=var)
  if not CheckHistoCategories(plt.hists[var], categories):
    print("Nope")
    return
  plt.ResetExtraBkg(0)
  plt.SetRatio(doRatio)
  plt.SetOutpath(outpath)
  plt.plotData = doData
  plt.SetLumi(lumi, "pb$^{-1}$", '13.6 TeV')
  plt.SetRatioRange(0.75,1.25)

  if var in ['counts', 'njets','nbtags', 'met', 'jet0pt', 'jet0eta', 'jetpt', 'jeteta', 'ht', 'invmass'] and not 'sign' in categories:
    categories['sign'] = 'OS'
  #plt.SetDataName('Pseudodata')
  if nonpromptDD:
    hNonprompt = GetNonprompt(plt, var, categories)
    plt.AddExtraBkgHist(hNonprompt)

  label = (GetChLab(categories['channel']) if isinstance(categories['channel'], str) else GetChLab(categories['channel'][0]) ) + GetLevLab(categories['level']) + (' SS' if categories['sign'] == 'SS' else '')
  plt.SetCategories(categories)
  #AddLabel(self, x, #y, text, options={}):
  plt.SetRegion(label)
  plt.SetOutput(output)

  b0 = None; bN = None
  if   var == 'deltaphi':
    b0 = 2
  elif var == 'invmass':
    b0 = 2

  if b0 is not None:
    plt.SetRebin(var, b0, bN, includeLower=True, includeUpper=True)
  
  plt.Stack(var, xtit='', ytit='', dosyst=True)


outpath = outpatho
outpath = '/nfs/fanae/user/juanr/www/public/ttrun3/ttLepEFF/28oct2022/' 

def Print2lplots(SS=False):
  for c in ['em']:
    for l in ['incl', 'g2jets', 'g1btag']:
      outp = outpath+'/'+l+'/' + ('SS/' if SS else '')
      cat = {'channel':c, 'level':l}
      if SS: cat['sign'] = 'SS'
      for var in ['counts', 'njets', 'met',  'jet0pt', 'jet0eta', 'jetpt', 'jeteta','ht','invmass']:
        outname = "%s_%s_%s"%(var, c, l) if isinstance(c, str) else "%s_%s_%s"%(var, c[0], l)
        Draw(var, cat, outname, outpath=outp, nonpromptDD= (var!='invmass' and not SS))

def PrintTPplots():
  global outpath 
  outpath += '/TP/'
  for c in ['eprobe', 'muprobe']:
    for l in ['incl', 'g2jets']:
      outp = outpath+'/'+l+'/'
      for sign in ['OS', 'SS']:
        for lep in ['pass', 'fail']:
          cat = {'channel':c, 'level':l, 'sign':sign, 'lepton':lep if lep == 'pass' else ['pass', 'fail']}
          for var in ['lpt', 'leta']:
            outname = "%s_%s_%s_%s_%s"%(var, c[0], l, sign, lep if lep == 'pass' else 'all')
            Draw(var, cat, outname, outpath=outp, nonpromptDD=(sign=='OS'))

if not var is None:
  categories = { 'channel': ch, 'level' : level}#, 'syst':'norm'}
  outp = outpath+'/'+level+'/'
  outname = "%s_%s_%s"%(var, ch, level) if isinstance(ch, str) else "%s_%s_%s"%(var, ch[0], level)
  Draw(var, categories, outname, outpath=outp, nonpromptDD=True)

else:
  #Draw('njets', {'channel':['em'], 'level':'dilep'}, output='njetsmalnrm', outpath=outpath)
  Print2lplots(1)
  #PrintTPplots()

