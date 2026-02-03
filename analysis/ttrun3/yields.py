from config import *

def Draw(var, categories, output=None, label='', outpath='temp/', doRatio=True):
  plt = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi, var=var)
  if not CheckHistoCategories(plt.hists[var], categories):
    print("Nope")
    return
  plt.SetRatio(doRatio)
  plt.SetOutpath(outpath)
  plt.plotData = doData
  plt.SetLumi(lumi, "pb$^{-1}$", '13.6 TeV')
  plt.SetCategories(categories)
  plt.SetDataName('Pseudodata')
  label = (GetChLab(categories['channel']) if isinstance(categories['channel'], str) else GetChLab(categories['channel'][0]) ) + GetLevLab(categories['level'])
  #AddLabel(self, x, #y, text, options={}):
  plt.SetRegion(label)
  plt.SetOutput(output)
  
  plt.SetSystematics(syst=['lepSF', 'JES', 'trigSF', 'FSR', 'ISR'])#, 'lepSF']) # FSR, ISR, JES, lepSF, trigSF
  plt.Stack(var, xtit='', ytit='', dosyst=True)

  #plt.PrintYields('counts')

def Print2lplots():
  for c in ['em', 'ee', 'mm']:
    for l in ['dilep', 'g2jets']:
      outp = outpath+'/'+l+'/'
      cat = {'channel':c, 'level':l}#, 'syst':'norm'}
      for var in ['njets', 'ht', 'met', 'j0pt', 'j0eta', 'invmass', 'invmass2']:
        if l=='dilep' and var in ['j0pt', 'j0eta']: continue
        outname = "%s_%s_%s"%(var, c, l)
        Draw(var, cat, outname, outpath=outp)
      for var in ['counts', 'l0pt','ept', 'mpt', 'l0eta', 'eeta', 'meta']:
        cat['sign'] = 'OS'
        outname = "%s_%s_%s"%(var, c, l)
        Draw(var, cat, outname, outpath=outp)

outpath = '/nfs/fanae/user/juanr/www/public/ttrun3/usingRun2/' + outpatho
if not var is None:
  categories = { 'channel': ch, 'level' : level}#, 'syst':'norm'}
  outp = outpath+'/'+level+'/'
  outname = "%s_%s_%s"%(var, *ch, level)
  Draw(var, categories, outname, outpath=outp)


else:
  #Draw('invmass', {'channel':['em'], 'level':'g2jets', 'syst':'norm'}, output='invmass_em_g2jets', outpath=outpath)
  Print2lplots()


