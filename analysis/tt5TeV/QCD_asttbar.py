'''
 Estimate the QCD background
 Usagel: 
   qcd = QCD(path, prDic=processDic, bkglist=bkglist, lumi=lumi)
   hqcd = qcd.GetQCD(var, categories)
'''

from analysis.tt5TeV.config import *

class QCD:

  def __init__(self, path, prDic, bkglist, lumi, categories={}, var=None):
    self.path = path
    self.prDic = prDic
    self.bkglist = bkglist
    self.lumi = lumi
    self.varlist = [var] if isinstance(var, str) else var

    self.categories = categories
    #if 'channel' in self.categories:
    #  self.chan = self.categories['channel']
    #  self.categories.pop('channel')
    self.cat_fake = self.categories.copy()
    self.doAllChan = False
    if not 'channel' in self.categories:
      #self.doAllChan = True
      #self.categories['channel'] = ['e', 'm']
      #self.cat_fake  ['channel'] = ['e_fake', 'm_fake']
      pass
    else:
      self.cat_fake['channel'] = [(c+'_fake' if not 'fake' in c else c) for c in self.cat_fake['channel']]

    # Initialize histograms and normalization 
    self.QCDhist = {}
    self.norm = {-1 : 1, 0:1, -1:1}
    self.LoadHistos()

  def LoadQCDForVar(self, var, group=False):
    ''' Compute (data - MC) for a given distribution (dense_axis) '''
    if self.doAllChan:
      h_data_fake = self.plt.GetHistogram(var, ['data'],{}).group('channel', hist.Cat("channel", "channel"), {'e':'e_fake', 'm':'m_fake'}).group('process', hist.Cat("process", "process"), {'QCD':'data'} )
      h_mc_fake   = self.plt.GetHistogram(var, self.bkglist, {}).group('channel', hist.Cat("channel", "channel"), {'e':'e_fake', 'm':'m_fake'}).group('process', hist.Cat("process", "process"), {'QCD':self.bkglist})
    else:
      if 'channel' in self.categories:
        cat_fake = self.categories.copy()
        cat_fake['channel'] = [(c+'_fake' if not 'fake' in c else c) for c in cat_fake['channel']]
      else:
        cat_fake = self.cat_fake
      h_data_fake = self.plt.GetHistogram(var, ['data'], self.cat_fake, keepCats=True).group('process', hist.Cat("process", "process"), {'QCD':'data'})
      h_mc_fake   = self.plt.GetHistogram(var, self.bkglist, self.cat_fake, keepCats=True).group('process', hist.Cat("process", "process"), {'QCD':self.bkglist})
    h_mc_fake.scale(-1*self.lumi)
    systlist = [x.name for x in list(h_mc_fake.identifiers('syst'))]
    h_data_fake = self.FillDataSystCategories(h_data_fake, systlist, var)
    htot = (h_data_fake + h_mc_fake)
    if group:
      for cat in self.categories:
        htot = htot.group(cat, hist.Cat(cat, cat), {self.categories[cat]:self.categories[cat]})
    return htot

  def LoadHistos(self):
    ''' Load all the histograms '''
    self.plt = plotter(self.path, prDic=self.prDic,  bkgList=self.bkglist, lumi=self.lumi)
    #self.cat_fake = self.categories.copy()
    if self.varlist is None: self.varlist = self.plt.GetListOfVars()
    for var in self.varlist:
      if not CheckHistoCategories(self.plt.GetHistogram(var), {'channel' : 'e_fake', 'process':['data']+self.bkglist}, checkValues=True) and not CheckHistoCategories(self.plt.GetHistogram(var), {'channel' : 'm_fake', 'process':['data']+self.bkglist}, checkValues=True):
        print("skipping var: ", var)
        continue
      dense_axes = [x.name for x in self.plt.GetHistogram(var).dense_axes()]
      if len(dense_axes) > 1: continue
      if not var == dense_axes[0]: continue
      if var in self.QCDhist.keys():
        self.QCDhist[var] += self.LoadQCDForVar(var)
      else:
        self.QCDhist[var] = self.LoadQCDForVar(var)

  def FillDataSystCategories(self, hdata, systlist, var):
    ''' Fill systematics for fake data distribution (using nominal) '''
    hnorm = hdata.integrate('process', 'QCD')
    if CheckHistoCategories(hnorm, {'syst':'norm'}):
      hnorm = hnorm.integrate('syst', 'norm')
    levels = [x.name for x in list(hdata.identifiers('level'))]
    channels = ['e', 'm']
    for l in levels:
      for c in channels:
        hnorm2 = hnorm.integrate('level', l).integrate('channel', c)
        if hnorm2.values() == {}: continue
        for s in systlist:
          if s == 'norm': continue
          bins, vals = GetXYfromH1D(hnorm2, axis=var, mode='centers', errors=False, overflow=False)
          hdata.fill(**{'syst':s, 'weight':vals, 'process':'QCD', 'channel':c, 'level':l, var:bins})
    '''
    return hdata
    if not 'level' in self.categories:
      levels = [x.name for x in list(hdata.identifiers('level'))]
      for l in levels:
        for c in ['e', 'm']:
          if not c in [x.name for x in hnorm.identifiers('channel')]: continue
          hnorm2 = hnorm.integrate('level', l).integrate('channel', c)
          if hnorm2.values() == {}: continue
          for s in systlist:
            if s == 'norm': continue
            bins, vals = GetXYfromH1D(hnorm2, axis=var, mode='centers', errors=False, overflow=False)
            hdata.fill(**{'syst':s, 'weight':vals, 'process':'QCD', 'channel':c, 'level':l, var:bins})
    else:
      print('hnorm')
      PrintHisto(hnorm)
      bins, vals = GetXYfromH1D(hnorm, axis=var, mode='centers', errors=False, overflow=False)
      for s in systlist:
        if s == 'norm': continue
        hdata.fill(**{'syst':s, 'weight':vals, 'process':'QCD', var:bins})
    '''
    return hdata

  def GetNormalization(self, categories, sys=0):
    ''' Load normalization '''
    chan = categories['channel']
    catfake = categories.copy()
    catfake['channel'] = [(c+'_fake' if not 'fake' in c else c) for c in catfake['channel']]
    if 'fake' in chan:
      categories['channel'] = chan.replace('_fake', '')
    else:
      catfake['channel'] = chan + '_fake'
    if   sys== 1: countsData = 'counts_metl25'
    elif sys==-1: countsData = 'counts_metl15'
    else:         countsData = 'counts_metl20'
    data_metl20    = self.plt.GetYields(countsData, categories, pr='data', overflow='none')
    mc_metl20      = self.plt.GetYields(countsData, categories, pr=self.bkglist, overflow='none')
    data_metl20_fk = self.plt.GetYields(countsData, catfake,   pr='data', overflow='none')
    mc_metl20_fk   = self.plt.GetYields(countsData, catfake,   pr=self.bkglist, overflow='none')
    fact = (data_metl20 - mc_metl20)/(data_metl20_fk - mc_metl20_fk)
    return fact

  def ScaleWithSyst(self, hqcd, categories, var):
    h   = self.QCDhist[var].copy()
    hup = h.copy()
    hdo = h.copy()
    fact   = self.GetNormalization(categories, 0)
    factup = self.GetNormalization(categories, 1)
    factdo = self.GetNormalization(categories,-1)
    h.scale(fact)
    hup.scale(factup)
    hdo.scale(factdo)
    hup = GroupKeepOrder(hup, [['syst', 'syst', {'QCDUp':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
    hdo = GroupKeepOrder(hdo, [['syst', 'syst', {'QCDDown':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
    h += (hup + hdo)
    return h

  def GetQCDall(self):
    ''' Get QCD histogram for all variables '''
    variables = list(self.QCDhist.keys())
    var = variables[0]
    h = self.QCDhist[var]
    ident = list(h.identifiers('level'))
    levels = [x.name for x in list(self.QCDhist[variables[0]].identifiers('level'))]
    channels = ['e', 'm']
    histos = {}
    total = len(variables)*len(channels)*len(levels)
    progress = 0
    for v in variables:
      for c in channels:
        for l in levels:
          # Cool progress bar every loop
          progress100 = (progress / total) * 100
          print("\r[{:<50}] {:.2f} % ({:.0f}/{:.0f} )".format('#' * int(progress100/2), progress100, progress, total), end="")
          self.categories = {'level':l, 'channel':c}
          h = self.LoadQCDForVar(v, group=True)
          self.ScaleWithSyst(h, self.categories, v)
          if v in list(histos.keys()): histos[v] += h
          else: histos[v] = h
          progress += 1
    return histos
    

  def GetNorm(self, sys=0):
    return self.norm[sys]

  def GetQCD(self, var, categories={}, sys=False):
    if categories=={}: categories = self.categories
    fact = self.GetNormalization(categories, 0)
    #print('sys = ', sys, 'fact = ', fact)
    h = self.QCDhist[var].copy()
    for cat in categories: 
      axes = list(h.sparse_axes())
      if not cat in axes: continue
      h = h.integrate(cat, categories[cat])
    if sys:
      hup = h.copy()
      hdo = h.copy()
      factup = self.GetNormalization(categories, 1)
      factdo = self.GetNormalization(categories,-1)
      hup.scale(factup)
      hdo.scale(factdo)
      hup = GroupKeepOrder(hup, [['syst', 'syst', {'QCDUp':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
      hdo = GroupKeepOrder(hdo, [['syst', 'syst', {'QCDDown':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
    h.scale(fact)
    if sys:
      h += (hup + hdo)
    #if sys==1: # Up
    #  #h = GroupKeepOrder(h, [['syst', 'syst', {'QCDUp':'norm'}, 'keep']])
    #  h = GroupKeepOrder(h, [['syst', 'syst', {'QCDUp':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
    #elif sys==-1: # Down
    #  #h = GroupKeepOrder(h, ['syst', 'syst', {'QCDDown':'norm'}])
    #  h = GroupKeepOrder(h, [['syst', 'syst', {'QCDDown':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
    #print("sys = ", sys, ", h = ", h)
    #PrintHisto(h)
    return h

  def GetYield(self, categories, var='counts'):
    hqcd = self.GetQCD(var, categories)
    for cat in categories:
      axes = list(hqcd.sparse_axes())
      if not cat in axes: continue
      hqcd = hqcd.integrate(cat, categories[cat])
    return hqcd.values()[('QCD',)][0]


if __name__ == '__main__':
  from config import *
  qcd = QCD(path, prDic=processDic_noQCD, bkglist=bkglist_noQCD, lumi=lumi, categories={'channel':ch, 'level':level}, var=var)
  hqcd = qcd.GetQCD(var, {'channel':ch, 'level':level}, 1)
  qcdnom = hqcd.integrate('syst', 'norm').integrate('process', 'QCD')
  qcdup  = hqcd.integrate('syst', 'QCDUp').integrate('process', 'QCD')
  qcddo  = hqcd.integrate('syst', 'QCDDown').integrate('process', 'QCD')

  labels = ['QCD', 'QCD Up', 'QCD Down']
  colors = ['black', 'red', 'blue']
  DrawComp([qcdnom, qcdup, qcddo], colors=colors, axis=var, title='', labels=labels, xtit=var, ytit=None, doFill=None, outname=output)
