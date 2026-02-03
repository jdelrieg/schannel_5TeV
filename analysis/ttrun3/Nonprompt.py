from config import *

prdic = {'Prompt': 'TTTo2L2Nu, tbarW, tW, WWTo2L2Nu, WZTo3LNu, DYJetsToLL_M50, DY_M10to50', 'Nonprompt': 'WJetsToLNu, TTToSemiLep'}

p = plotter(path, prDic=prdic, bkgList=bkglist, colors=colordic, lumi=lumi)
p.SetDataName('pseudodata')

def GetYieldsNP(plotter, sign='OS', process='Nonprompt', level='g2jets'):
  cat = {'level':level, 'syst':'norm', 'sign':sign, 'channel':'em'}
  if process == 'data':
    h = plotter.GetHistogram('counts', process=None, categories=cat)
    h.scale(lumi)
    h,_ = plotter.GetData('counts', h)
  else:
    h = plotter.GetHistogram('counts', process=process, categories=cat)
    h.scale(lumi)
  y, e = h.values(sumw2=True)[()]
  if process == 'data': e = np.sqrt(y)
  return y, e

def Nonprompt(plt, level='g2jets', save=False):
  NdataSS, NdataSS_e = GetYieldsNP(p, 'SS', 'data', level)
  NpSS   , NpSS_e    = GetYieldsNP(p, 'SS', 'Prompt', level)

  NnpOS, NnpOS_e = GetYieldsNP(p, 'OS', 'Nonprompt', level)
  NnpSS, NnpSS_e = GetYieldsNP(p, 'SS', 'Nonprompt', level)

  NSS, NSS_e = (NdataSS - NpSS, NdataSS_e - NpSS_e)
  fact   = NnpOS / NnpSS
  fact_e = fact*np.sqrt( (NnpOS_e/NnpOS)**2 + (NnpSS_e/NnpSS)**2 )

  NOS   = NSS*fact
  NOS_e = NOS * np.sqrt( (NSS_e/NSS)**2 + (fact_e/fact)**2 )

  if save:
    t = OutText(outpatho, 'Nonprompt_'+level, 'new', 'tex')
    t.bar()
    t.SetTexAlign('l c')
    t.line("Source" + t.vsep() + "$\mathrm{e}^\pm\mu^\mp$"); t.sep()
    t.line("Prompt SS (MC)" + t.vsep() + "%1.2f $\pm$ %1.2f"%(NpSS, NpSS_e));
    t.line("Data SS" + t.vsep() + "%1.2f $\pm$ %1.2f"%(NdataSS, NdataSS_e));
    t.line("Data - prompt (MC) SS" + t.vsep() + "%1.2f $\pm$ %1.2f"%(NSS, NSS_e)); t.sep()
    
    t.line("Nonprompt SS (MC)" + t.vsep() + "%1.2f $\pm$ %1.2f"%(NnpSS, NnpSS_e) );
    t.line("Nonprompt OS (MC)" + t.vsep() + "%1.2f $\pm$ %1.2f"%(NnpOS, NnpOS_e) );
    t.line("Ratio OS/SS (MC)" + t.vsep() + "%1.2f $\pm$ %1.2f"%(fact, fact_e) ); t.sep()

    t.line("Nonprompt estimate" + t.vsep() + "%1.2f $\pm$ %1.2f"%(NOS, NOS_e) ); t.sep()
    t.write()

  return NOS, NOS_e

def DrawNonprompt(p, var='l0pt', level=level, chan=ch):
  #prdic =  #{'tt': 'TTTo2L2Nu, tbarW, tW', 'DY': 'WWTo2L2Nu, WZTo3LNu, DYJetsToLL_M50, DY_M10to50', 'Nonprompt': 'WJetsToLNu, TTToSemiLep'}
  #colordic = {'tt' : '#ff6782', 'DY':'#00b4bd', 'Nonprompt': '#47ce33'}
  prdic = {'Prompt': 'TTTo2L2Nu, tbarW, tW, WWTo2L2Nu, WZTo3LNu, DYJetsToLL_M50, DY_M10to50', 'Nonprompt': 'WJetsToLNu, TTToSemiLep'}
  colordic = {'Prompt' : '#ff6782', 'Nonprompt': '#47ce33'}
  bkgList = ['Prompt','Nonprompt']
  p = plotter(path, prDic=prdic, bkgList=bkgList, colors=colordic, lumi=lumi, var=var)
  p.SetDataName('pseudodata')
  categories = {'level':level, 'sign':'SS', 'channel':chan, 'syst':'norm'}
  p.SetRatio(True)
  p.SetOutpath(outpatho)
  p.plotData = True
  p.SetLumi(lumi, "pb$^{-1}$", '13.6 TeV')
  p.SetCategories(categories)
  label = (GetChLab(categories['channel']) if isinstance(categories['channel'], str) else GetChLab(categories['channel'][0]) ) + GetLevLab(categories['level'])  + ', same-sign'
  p.SetRebin(var, 2)
  p.SetRegion(label)
  p.SetOutput(output)
  #p.SetSystematics(syst=['lepSF', 'JES', 'trigSF', 'FSR', 'ISR'])#, 'lepSF']) # FSR, ISR, JES, lepSF, trigSF
  p.Stack(var, xtit='', ytit='', dosyst=True)



if __name__ == '__main__':
  #Nonprompt(p, level=level, save=True)
  DrawNonprompt(p, 'njets')


