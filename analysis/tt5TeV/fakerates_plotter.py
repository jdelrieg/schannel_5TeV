'''
 Estimate the QCD background
 Usagel: 
   qcd = QCD(path, prDic=processDic, bkglist=bkglist, lumi=lumi)
   hqcd = qcd.GetQCD(var, categories)
'''

from analysis.tt5TeV.config import *
import time
import matplotlib.pyplot as plot


if os.path.exists(path + 'QCD.pkl.gz'):
  print('WARNING: QCD file already exists in path... moving to ', path + 'old/QCD.pkl.gz.old')
  if not os.path.exists(path + 'old/'):
    os.makedirs(path + 'old/')
  os.rename(path + 'QCD.pkl.gz', path + 'old/QCD.pkl.gz.old')

print('Loading files from ', path, '... (this may take a while)') 
plt = plotter(path, prDic=processDic_noQCD,  bkgList=bkglist_noQCD, lumi=lumi)

# Get list of variables
varlist = plt.GetListOfVars()
variables = []
print('Getting list of variables for which QCD can be estimated...')
print("  -- skipping variables: ", var, end='')
for var in varlist:
  if not CheckHistoCategories(plt.GetHistogram(var), {'channel' : 'e_fake_plus', 'process':['data']+bkglist_noQCD}, checkValues=True) and not CheckHistoCategories(plt.GetHistogram(var), {'channel' : 'm_fake_plus', 'process':['data']+bkglist_noQCD}, checkValues=True) and not CheckHistoCategories(plt.GetHistogram(var), {'channel' : 'm_fake_minus', 'process':['data']+bkglist_noQCD}, checkValues=True) and not CheckHistoCategories(plt.GetHistogram(var), {'channel' : 'e_fake_minus', 'process':['data']+bkglist_noQCD}, checkValues=True):
    print(', ', var, end='')
    continue
  variables.append(var)

# Channels and levels
print('Getting levels and systematics...')
channels    = ['e_plus','e_minus', 'm_plus','m_minus']
#channels    = ['m_plus','m_minus']
levels      = [x.name for x in list(plt.GetHistogram(variables[0]).identifiers('level'))]
systematics = [x.name for x in list(plt.GetHistogram(variables[0]).identifiers('syst'))]
if 'norm' in systematics: systematics.remove('norm')
systematics=[]

print('--------------------------------------------------')
print('Variables:', variables)
print('Channels:', channels)
print('Levels:', levels)
print('Systematics:', systematics)
print('--------------------------------------------------')

def GetQCDforVar(var):
  ''' Compute (data - MC) for a given variable (dense_axis) -- all channels and levels!! '''
  catfake = {'channel' : ['e_fake_plus','e_fake_minus', 'm_fake_plus','m_fake_minus']}
  h_data_fake = plt.GetHistogram(var, ['data']     , catfake, keepCats=True)#.group('channel', hist.Cat("channel", "channel"), {'e_fake':'e', 'm_fake':'m'}).group('process', hist.Cat("process", "process"), {'QCD':'data'} )
  h_mc_fake   = plt.GetHistogram(var, bkglist_noQCD, catfake, keepCats=True)#.group('channel', hist.Cat("channel", "channel"), {'e_fake':'e', 'm_fake':'m'}).group('process', hist.Cat("process", "process"), {'QCD':bkglist_noQCD})
  h_data_fake = GroupKeepOrder(h_data_fake, [['channel', 'channel', {'e_plus':'e_fake_plus','e_minus':'e_fake_minus', 'm_plus':'m_fake_plus','m_minus':'m_fake_minus'}], ['process', 'process', {'QCD':'data'       }]])
  h_mc_fake   = GroupKeepOrder(h_mc_fake  , [['channel', 'channel', {'e_plus':'e_fake_plus','e_minus':'e_fake_minus', 'm_plus':'m_fake_plus','m_minus':'m_fake_minus'}], ['process', 'process', {'QCD':bkglist_noQCD}]])
  
  h_mc_fake.scale(-1*lumi)
  FillDataSystCategories(h_data_fake)#, var=var)
  htot = (h_data_fake + h_mc_fake)
  hd  = h_data_fake.integrate('channel', 'e_plus').integrate('process').integrate('level', '2j1b').integrate('syst', 'norm') #4j2b
  hmc = h_mc_fake.integrate('channel', 'e_plus').integrate('process').integrate('level', '2j1b').integrate('syst', 'norm')   #4j2b
  ht  = htot.integrate('channel', 'e_plus').integrate('process').integrate('level', '2j1b').integrate('syst', 'norm')        #4j2b
  if htot._sumw != {}: 
    for i in htot._sumw.keys():
      htot._sumw[i]=np.where(htot._sumw[i]<0,0.0,htot._sumw[i])
  return htot

def GroupCats(h, catdic):
  ''' Group categories '''
  for cat in catdic:
    h = h.group(cat, hist.Cat(cat, cat), {catdic[cat]:catdic[cat]})
  return h

def FillDataSystCategories(hdata, var=None):
  ''' Fill systematics for fake data distribution (using nominal) '''
  if var is None: var = hdata.dense_axes()[0].name
  hnorm = hdata.copy()
  if CheckHistoCategories(hnorm, {'syst':'norm'}):
    hnorm = hnorm.integrate('syst', 'norm')
  if CheckHistoCategories(hnorm, {'process':'QCD'}):
    hnorm = hnorm.integrate('process', 'QCD')
  elif CheckHistoCategories(hnorm, {'process':'data'}):
    hnorm = hnorm.integrate('process', 'data')
  for l in levels:
    for c in channels:
      hnorm2 = hnorm.integrate('level', l).integrate('channel', c)
      if hnorm2.values() == {}: continue
      bins, vals = GetXYfromH1D(hnorm2, axis=var, mode='centers', errors=False, overflow=False)
      for s in systematics:
        if s == 'norm': continue
        hdata.fill(**{'syst':s, 'weight':vals, 'process':'QCD', 'channel':c, 'level':l, var:bins})

def GetQCDnorm(chan, level, sys=0):
  cat     = {'channel':chan, 'level':level}
  catfake = {'channel':chan[0] + '_fake_' + chan[2:], 'level':level}
  
  if   sys== 1: countsData = 'counts_metl30' #l25 original #l180  ht  #l30 mtw que uso en mtw30   #l10  met con corte en 15   45 
  elif sys==-1: countsData = 'counts_metl20' #l15			#l170      #l20						  #l5						  35
  else:         countsData = 'counts_metl25' #l20          #l175      #l25 						  #l75(7.5)					  40
  
  data_metl20    = plt.GetYields(countsData, cat    , pr='data'       , overflow='none')
  mc_metl20      = plt.GetYields(countsData, cat    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk = plt.GetYields(countsData, catfake, pr='data'       , overflow='none')
  mc_metl20_fk   = plt.GetYields(countsData, catfake, pr=bkglist_noQCD, overflow='none')
  
  #if sys==0:print(chan,level);print('data_metl20',data_metl20,'mc_metl20',mc_metl20,'data_metl20_fk',data_metl20_fk,'mc_metl20_fk',mc_metl20_fk)
  fact = (data_metl20 - mc_metl20)/(data_metl20_fk - mc_metl20_fk)
  #if sys==0:print('fact',chan,level,fact)
  return fact


def fakerates(chan, level):
  cat     = {'channel':chan, 'level':level}
  catfake = {'channel':chan[0] + '_fake_' + chan[2:], 'level':level}

  
  countsData = 'counts_metl25_loweta_lowpt' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_metl20    = plt.GetYields(countsData, cat    , pr='data'       , overflow='none')
  mc_metl20      = plt.GetYields(countsData, cat    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk = plt.GetYields(countsData, catfake, pr='data'       , overflow='none')
  mc_metl20_fk   = plt.GetYields(countsData, catfake, pr=bkglist_noQCD, overflow='none')
  factll = (data_metl20 - mc_metl20)/(data_metl20_fk - mc_metl20_fk)
  print(chan,level,'factll',factll,'a',data_metl20,'b',mc_metl20,'c',data_metl20_fk,'d',mc_metl20_fk)
  
  countsData = 'counts_metl25_loweta_highpt' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_metl20    = plt.GetYields(countsData, cat    , pr='data'       , overflow='none')
  mc_metl20      = plt.GetYields(countsData, cat    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk = plt.GetYields(countsData, catfake, pr='data'       , overflow='none')
  mc_metl20_fk   = plt.GetYields(countsData, catfake, pr=bkglist_noQCD, overflow='none')
  factlh = (data_metl20 - mc_metl20)/(data_metl20_fk - mc_metl20_fk)
  print(chan,level,'factlh',factlh,'a',data_metl20,'b',mc_metl20,'c',data_metl20_fk,'d',mc_metl20_fk)

  countsData = 'counts_metl25_higheta_lowpt' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_metl20    = plt.GetYields(countsData, cat    , pr='data'       , overflow='none')
  mc_metl20      = plt.GetYields(countsData, cat    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk = plt.GetYields(countsData, catfake, pr='data'       , overflow='none')
  mc_metl20_fk   = plt.GetYields(countsData, catfake, pr=bkglist_noQCD, overflow='none')
  facthl = (data_metl20 - mc_metl20)/(data_metl20_fk - mc_metl20_fk)
  print(chan,level,'facthl',facthl,'a',data_metl20,'b',mc_metl20,'c',data_metl20_fk,'d',mc_metl20_fk)

  countsData = 'counts_metl25_higheta_highpt' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_metl20    = plt.GetYields(countsData, cat    , pr='data'       , overflow='none')
  mc_metl20      = plt.GetYields(countsData, cat    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk = plt.GetYields(countsData, catfake, pr='data'       , overflow='none')
  mc_metl20_fk   = plt.GetYields(countsData, catfake, pr=bkglist_noQCD, overflow='none')
  facthh = (data_metl20 - mc_metl20)/(data_metl20_fk - mc_metl20_fk)
  print(chan,level,'facthh',facthh,'a',data_metl20,'b',mc_metl20,'c',data_metl20_fk,'d',mc_metl20_fk)

  return factll,factlh,facthl,facthh

def GetQCDnorm_v2(chan, level, sys=0):
  cat     = {'channel':chan, 'level':level}
  catfake = {'channel':chan[0] + '_fake_' + chan[2:], 'level':level}
  countsData = 'counts' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_fk = plt.GetYields(countsData, catfake, pr='data'       , overflow='none')
  mc_fk   = plt.GetYields(countsData, catfake, pr=bkglist_noQCD, overflow='none')
  
  if sys==0:print(chan,level);print('data_fk',data_fk,'mc_fk',mc_fk)
  fact = (data_fk - mc_fk)
  if sys==0:print('fact old',chan,level,fact)
  return fact


def NormQCD(hqcd, chan, level):
  ''' Normalize QCD '''
  fact   = GetQCDnorm(chan, level)
  #otro_factor=GetQCDnorm_v2(chan,level)
  factUp = GetQCDnorm(chan, level, sys=1)
  factDo = GetQCDnorm(chan, level, sys=-1)
  if fact<0:
    fact=0.0; factUp=0.0; factDo=0.0
  cat = {'channel':chan, 'level':level}
  #for c in cat:
  #  hqcd = hqcd.integrate(c, cat[c])
  hqcd = GroupKeepOrder(hqcd, [['channel', 'channel', {cat['channel']:cat['channel']}], ['level', 'level', {cat['level']:cat['level']}]])
  hqcdUp = hqcd.copy()
  hqcdDo = hqcd.copy()
  #print(hqcd.sum(),hqcd.values())
  hqcd  .scale(fact)
  hqcdUp.scale(factUp)
  hqcdDo.scale(factDo)
  hqcdUp = GroupKeepOrder(hqcdUp, [['syst', 'syst', {'QCDUp':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
  hqcdDo = GroupKeepOrder(hqcdDo, [['syst', 'syst', {'QCDDown':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
  hqcd += hqcdUp
  hqcd += hqcdDo
  return hqcd

#def GetQCD(qcdHist, level, chan):
#  ''' Apply QCD normalization safely '''
#  hqcd = qcdHist.copy()
#  #hqcd = h.group('level', hist.Cat("level", "level"), {level:level}).group('channel', hist.Cat("channel", "channel"), {chan:chan})
  #GroupKeepOrder(hqcd, [['level', 'level', {level:level}], ['channel', 'channel', {chan:chan}] ])
#  hqcd = NormQCD(hqcd, chan, level)
#  return hqcd

def GetQCDpar(inputs):
  ''' To be used with multiprocessing '''
  qcdHist, var, level, chan, outdict = inputs
  h = NormQCD(qcdHist, chan, level)
  h = GroupKeepOrder(h, [['process', 'sample', {'QCD':'QCD'}]])
  if var not in outdict: outdict[var] = h
  else: outdict[var] += h
  outdict['progress'] += 1
  total = outdict['total']
  print("\r[{:<50}] {:.2f} %".format('#' * int(outdict['progress']/total*100/2), float(outdict['progress'])/total*100), end='')
  
################################################ RUN
# Initialize histograms
histos = {}
d0 = time.time()
progress = 0
total = len(variables)*len(levels)*len(channels)
QCDhistos = {}
print('Calculating distributions...')
ivar = 0; 
for var in variables:
  progress100 = (ivar / len(variables)) * 100
  print("\r[{:<50}] {:.2f} %".format('#' * int(progress100/2), progress100), end='')
  QCDhistos[var] = GetQCDforVar(var)
  # Normalize back to 1/pb
  QCDhistos[var].scale(1./lumi)
  ivar += 1
print("\r[{:<50}] {:.2f} %".format('#' * int(progress100/2), progress100))

print('Grouping and normalizing...')
from multiprocessing import Pool, Manager
pool = Pool(nSlots)
manager = Manager()
outdict = manager.dict()
outdict['progress'] = 0
outdict['total'] = total

print('out',outdict)
# Create inputs
inputs = []
for var in variables:
  for l in levels:
    for c in channels:
      inputs.append([QCDhistos[var], var, l, c, outdict])

# Run in parallel
print('Running in parallel with {} slots...'.format(nSlots))
pool.map(GetQCDpar, inputs)
pool.close()
pool.join()

# Save histograms
#outdict = dict(outdict)
#saveHistos(path, 'QCD', outdict)

#print('\nQCD saved to ', path+'/QCD.pkl.gz') 

def fakerates_incl(chan):
  cat     = {'channel':chan+'_minus', 'level':'incl_tch'}
  catfake = {'channel':chan+'_fake_minus', 'level':'incl_tch'}

  cat1     = {'channel':chan+'_plus', 'level':'incl_tch'}
  catfake1 = {'channel':chan+'_fake_plus', 'level':'incl_tch'}

  countsData = 'counts_metl25_loweta_lowpt' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_metl20    = plt.GetYields(countsData, cat    , pr='data'       , overflow='none')
  mc_metl20      = plt.GetYields(countsData, cat    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk = plt.GetYields(countsData, catfake, pr='data'       , overflow='none')
  mc_metl20_fk   = plt.GetYields(countsData, catfake, pr=bkglist_noQCD, overflow='none')
  
  data_metl201    = plt.GetYields(countsData, cat1    , pr='data'       , overflow='none')
  mc_metl201      = plt.GetYields(countsData, cat1    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk1 = plt.GetYields(countsData, catfake1, pr='data'       , overflow='none')
  mc_metl20_fk1   = plt.GetYields(countsData, catfake1, pr=bkglist_noQCD, overflow='none')
  factll = (data_metl20+data_metl201 - mc_metl20-mc_metl201)/(data_metl20_fk+data_metl20_fk1 - mc_metl20_fk-mc_metl20_fk1)
  print('factll',factll,'a',data_metl20+data_metl201,'b',mc_metl20+mc_metl201,'c',data_metl20_fk+data_metl20_fk1,'d',mc_metl20_fk+mc_metl20_fk1)
  
  countsData = 'counts_metl25_loweta_highpt' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_metl20    = plt.GetYields(countsData, cat    , pr='data'       , overflow='none')
  mc_metl20      = plt.GetYields(countsData, cat    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk = plt.GetYields(countsData, catfake, pr='data'       , overflow='none')
  mc_metl20_fk   = plt.GetYields(countsData, catfake, pr=bkglist_noQCD, overflow='none')
  data_metl201    = plt.GetYields(countsData, cat1    , pr='data'       , overflow='none')
  mc_metl201      = plt.GetYields(countsData, cat1    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk1 = plt.GetYields(countsData, catfake1, pr='data'       , overflow='none')
  mc_metl20_fk1   = plt.GetYields(countsData, catfake1, pr=bkglist_noQCD, overflow='none')  
  factlh = (data_metl20+data_metl201 - mc_metl20-mc_metl201)/(data_metl20_fk+data_metl20_fk1 - mc_metl20_fk-mc_metl20_fk1)
  print('factlh',factlh,'a',data_metl20+data_metl201,'b',mc_metl20+mc_metl201,'c',data_metl20_fk+data_metl20_fk1,'d',mc_metl20_fk+mc_metl20_fk1)	
 	

  countsData = 'counts_metl25_higheta_lowpt' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_metl20    = plt.GetYields(countsData, cat    , pr='data'       , overflow='none')
  mc_metl20      = plt.GetYields(countsData, cat    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk = plt.GetYields(countsData, catfake, pr='data'       , overflow='none')
  mc_metl20_fk   = plt.GetYields(countsData, catfake, pr=bkglist_noQCD, overflow='none')
  data_metl201    = plt.GetYields(countsData, cat1    , pr='data'       , overflow='none')
  mc_metl201      = plt.GetYields(countsData, cat1    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk1 = plt.GetYields(countsData, catfake1, pr='data'       , overflow='none')
  mc_metl20_fk1   = plt.GetYields(countsData, catfake1, pr=bkglist_noQCD, overflow='none')  
  facthl = (data_metl20+data_metl201 - mc_metl20-mc_metl201)/(data_metl20_fk+data_metl20_fk1 - mc_metl20_fk-mc_metl20_fk1)
  print('facthl',facthl,'a',data_metl20+data_metl201,'b',mc_metl20+mc_metl201,'c',data_metl20_fk+data_metl20_fk1,'d',mc_metl20_fk+mc_metl20_fk1)

  countsData = 'counts_metl25_higheta_highpt' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_metl20    = plt.GetYields(countsData, cat    , pr='data'       , overflow='none')
  mc_metl20      = plt.GetYields(countsData, cat    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk = plt.GetYields(countsData, catfake, pr='data'       , overflow='none')
  mc_metl20_fk   = plt.GetYields(countsData, catfake, pr=bkglist_noQCD, overflow='none')
  data_metl201    = plt.GetYields(countsData, cat1    , pr='data'       , overflow='none')
  mc_metl201      = plt.GetYields(countsData, cat1    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk1 = plt.GetYields(countsData, catfake1, pr='data'       , overflow='none')
  mc_metl20_fk1   = plt.GetYields(countsData, catfake1, pr=bkglist_noQCD, overflow='none')  
  facthh = (data_metl20+data_metl201 - mc_metl20-mc_metl201)/(data_metl20_fk+data_metl20_fk1 - mc_metl20_fk-mc_metl20_fk1)
  print('facthh',facthh,'a',data_metl20+data_metl201,'b',mc_metl20+mc_metl201,'c',data_metl20_fk+data_metl20_fk1,'d',mc_metl20_fk+mc_metl20_fk1)

  return factll,factlh,facthl,facthh

def fakerates_incl_total(chan):
  cat_e     = {'channel':'e_minus', 'level':'incl_tch'}
  catfake_e = {'channel':'e_fake_minus', 'level':'incl_tch'}

  cat_e1     = {'channel':'e_plus', 'level':'incl_tch'}
  catfake_e1 = {'channel':'e_fake_plus', 'level':'incl_tch'}

  cat_m     = {'channel':'m_minus', 'level':'incl_tch'}
  catfake_m = {'channel':'m_fake_minus', 'level':'incl_tch'}

  cat_m1     = {'channel':'m_plus', 'level':'incl_tch'}
  catfake_m1 = {'channel':'m_fake_plus', 'level':'incl_tch'}

  countsData = 'counts_metl25_loweta_lowpt' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_metl20    = plt.GetYields(countsData, cat_e    , pr='data'       , overflow='none')
  mc_metl20      = plt.GetYields(countsData, cat_e    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk = plt.GetYields(countsData, catfake_e, pr='data'       , overflow='none')
  mc_metl20_fk   = plt.GetYields(countsData, catfake_e, pr=bkglist_noQCD, overflow='none')
  
  data_metl201    = plt.GetYields(countsData, cat_e1    , pr='data'       , overflow='none')
  mc_metl201      = plt.GetYields(countsData, cat_e1    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk1 = plt.GetYields(countsData, catfake_e1, pr='data'       , overflow='none')
  mc_metl20_fk1   = plt.GetYields(countsData, catfake_e1, pr=bkglist_noQCD, overflow='none')

  data_metl20_m    = plt.GetYields(countsData, cat_m    , pr='data'       , overflow='none')
  mc_metl20_m      = plt.GetYields(countsData, cat_m    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk_m = plt.GetYields(countsData, catfake_m, pr='data'       , overflow='none')
  mc_metl20_fk_m   = plt.GetYields(countsData, catfake_m, pr=bkglist_noQCD, overflow='none')
  
  data_metl201_m    = plt.GetYields(countsData, cat_m1    , pr='data'       , overflow='none')
  mc_metl201_m      = plt.GetYields(countsData, cat_m1    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk1_m = plt.GetYields(countsData, catfake_m1, pr='data'       , overflow='none')
  mc_metl20_fk1_m   = plt.GetYields(countsData, catfake_m1, pr=bkglist_noQCD, overflow='none')
  factll = (data_metl20+data_metl201+data_metl20_m+data_metl201_m - mc_metl20-mc_metl201- mc_metl20_m-mc_metl201_m)/(data_metl20_fk+data_metl20_fk1+data_metl20_fk_m+data_metl20_fk1_m - mc_metl20_fk-mc_metl20_fk1- mc_metl20_fk_m-mc_metl20_fk1_m)
 
  countsData = 'counts_metl25_loweta_highpt' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_metl20    = plt.GetYields(countsData, cat_e    , pr='data'       , overflow='none')
  mc_metl20      = plt.GetYields(countsData, cat_e    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk = plt.GetYields(countsData, catfake_e, pr='data'       , overflow='none')
  mc_metl20_fk   = plt.GetYields(countsData, catfake_e, pr=bkglist_noQCD, overflow='none')
  
  data_metl201    = plt.GetYields(countsData, cat_e1    , pr='data'       , overflow='none')
  mc_metl201      = plt.GetYields(countsData, cat_e1    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk1 = plt.GetYields(countsData, catfake_e1, pr='data'       , overflow='none')
  mc_metl20_fk1   = plt.GetYields(countsData, catfake_e1, pr=bkglist_noQCD, overflow='none')

  data_metl20_m    = plt.GetYields(countsData, cat_m    , pr='data'       , overflow='none')
  mc_metl20_m      = plt.GetYields(countsData, cat_m    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk_m = plt.GetYields(countsData, catfake_m, pr='data'       , overflow='none')
  mc_metl20_fk_m   = plt.GetYields(countsData, catfake_m, pr=bkglist_noQCD, overflow='none')
  
  data_metl201_m    = plt.GetYields(countsData, cat_m1    , pr='data'       , overflow='none')
  mc_metl201_m      = plt.GetYields(countsData, cat_m1    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk1_m = plt.GetYields(countsData, catfake_m1, pr='data'       , overflow='none')
  mc_metl20_fk1_m   = plt.GetYields(countsData, catfake_m1, pr=bkglist_noQCD, overflow='none')
  factlh = (data_metl20+data_metl201+data_metl20_m+data_metl201_m - mc_metl20-mc_metl201- mc_metl20_m-mc_metl201_m)/(data_metl20_fk+data_metl20_fk1+data_metl20_fk_m+data_metl20_fk1_m - mc_metl20_fk-mc_metl20_fk1- mc_metl20_fk_m-mc_metl20_fk1_m)

  
  countsData = 'counts_metl25_higheta_lowpt' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_metl20    = plt.GetYields(countsData, cat_e    , pr='data'       , overflow='none')
  mc_metl20      = plt.GetYields(countsData, cat_e    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk = plt.GetYields(countsData, catfake_e, pr='data'       , overflow='none')
  mc_metl20_fk   = plt.GetYields(countsData, catfake_e, pr=bkglist_noQCD, overflow='none')
  
  data_metl201    = plt.GetYields(countsData, cat_e1    , pr='data'       , overflow='none')
  mc_metl201      = plt.GetYields(countsData, cat_e1    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk1 = plt.GetYields(countsData, catfake_e1, pr='data'       , overflow='none')
  mc_metl20_fk1   = plt.GetYields(countsData, catfake_e1, pr=bkglist_noQCD, overflow='none')

  data_metl20_m    = plt.GetYields(countsData, cat_m    , pr='data'       , overflow='none')
  mc_metl20_m      = plt.GetYields(countsData, cat_m    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk_m = plt.GetYields(countsData, catfake_m, pr='data'       , overflow='none')
  mc_metl20_fk_m   = plt.GetYields(countsData, catfake_m, pr=bkglist_noQCD, overflow='none')
  
  data_metl201_m    = plt.GetYields(countsData, cat_m1    , pr='data'       , overflow='none')
  mc_metl201_m      = plt.GetYields(countsData, cat_m1    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk1_m = plt.GetYields(countsData, catfake_m1, pr='data'       , overflow='none')
  mc_metl20_fk1_m   = plt.GetYields(countsData, catfake_m1, pr=bkglist_noQCD, overflow='none')
  facthl = (data_metl20+data_metl201+data_metl20_m+data_metl201_m - mc_metl20-mc_metl201- mc_metl20_m-mc_metl201_m)/(data_metl20_fk+data_metl20_fk1+data_metl20_fk_m+data_metl20_fk1_m - mc_metl20_fk-mc_metl20_fk1- mc_metl20_fk_m-mc_metl20_fk1_m)

  
  countsData = 'counts_metl25_higheta_highpt' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_metl20    = plt.GetYields(countsData, cat_e    , pr='data'       , overflow='none')
  mc_metl20      = plt.GetYields(countsData, cat_e    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk = plt.GetYields(countsData, catfake_e, pr='data'       , overflow='none')
  mc_metl20_fk   = plt.GetYields(countsData, catfake_e, pr=bkglist_noQCD, overflow='none')
  
  data_metl201    = plt.GetYields(countsData, cat_e1    , pr='data'       , overflow='none')
  mc_metl201      = plt.GetYields(countsData, cat_e1    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk1 = plt.GetYields(countsData, catfake_e1, pr='data'       , overflow='none')
  mc_metl20_fk1   = plt.GetYields(countsData, catfake_e1, pr=bkglist_noQCD, overflow='none')

  data_metl20_m    = plt.GetYields(countsData, cat_m    , pr='data'       , overflow='none')
  mc_metl20_m      = plt.GetYields(countsData, cat_m    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk_m = plt.GetYields(countsData, catfake_m, pr='data'       , overflow='none')
  mc_metl20_fk_m   = plt.GetYields(countsData, catfake_m, pr=bkglist_noQCD, overflow='none')
  
  data_metl201_m    = plt.GetYields(countsData, cat_m1    , pr='data'       , overflow='none')
  mc_metl201_m      = plt.GetYields(countsData, cat_m1    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk1_m = plt.GetYields(countsData, catfake_m1, pr='data'       , overflow='none')
  mc_metl20_fk1_m   = plt.GetYields(countsData, catfake_m1, pr=bkglist_noQCD, overflow='none')
  facthh = (data_metl20+data_metl201+data_metl20_m+data_metl201_m - mc_metl20-mc_metl201- mc_metl20_m-mc_metl201_m)/(data_metl20_fk+data_metl20_fk1+data_metl20_fk_m+data_metl20_fk1_m - mc_metl20_fk-mc_metl20_fk1- mc_metl20_fk_m-mc_metl20_fk1_m)

 
  return factll,factlh,facthl,facthh


factll,factlh,facthl,facthh=fakerates_incl('e')
print('inclusivos para electrones',factll,factlh,facthl,facthh)
factll,factlh,facthl,facthh=fakerates_incl('m')
print('inclusivos para muones',factll,factlh,facthl,facthh)
exit()


factors = {}

eta_values=np.array([-2.5, -0.8, 0.8, 2.5])
pt_values=np.array([0,80,200])


for channel in channels:
    factors[channel] = {}
    for level in levels:
        # Call your function to compute factors
        factll,factlh,facthl,facthh = fakerates(channel, level)#order is eta,pt
        

        # Assign factors to the dictionary
        factors[channel][level] = {
            "eta": {
                "values": eta_values,
                "factor_low_pt": np.array([facthl, factll, facthl]),#np.piecewise(eta_values,
                    #[eta_values < -0.8, (-0.8 <= eta_values) & (eta_values <= 0.8), eta_values > 0.8],
                    #[facthl, factll, facthl]
                #),
                "factor_high_pt": np.array([facthh, factlh, facthh]),#np.piecewise(eta_values,
                    #[eta_values < -0.8, (-0.8 <= eta_values) & (eta_values <= 0.8), eta_values > 0.8],
                    
                #),
            },
            "pt": {
                "values": pt_values,
                "factor_low_eta": np.array([factll, factlh]),#np.piecewise(pt_values,
                    #[pt_values < 90, pt_values >= 90],
                    
                #),
                "factor_high_eta": np.array([facthl, facthh]),#np.piecewise(pt_values,
                    #[pt_values < 90, pt_values >= 90],
                    #[facthl, facthh]
                #),
            }
        }


# Plotting function
def plot_factors(channel, level):
    if channel not in factors or level not in factors[channel]:
        print(f"Invalid channel '{channel}' or level '{level}'.")
        return

    level_data = factors[channel][level]

    # Plot for eta
    eta_data = level_data["eta"]
    plot.figure(figsize=(10, 6))
    #plot.plot(eta_data["values"], eta_data["factor_low_pt"], label="Low Pt", marker="o")
    #plot.plot(eta_data["values"], eta_data["factor_high_pt"], label="High Pt", marker="s")
    
		    
    for i in range(len(eta_data["values"])-1):  
        plot.hlines(eta_data["factor_low_pt"][i], eta_data["values"][i], eta_data["values"][i+1], colors='blue', linewidth=2,label=r'p$_{\mathrm{T}}<80$ GeV' if i == 0 else "")
    for i in range(len(eta_data["values"])-1):  
        plot.hlines(eta_data["factor_high_pt"][i], eta_data["values"][i], eta_data["values"][i+1], colors='red', linewidth=2,label=r'p$_{\mathrm{T}}\geq80$ GeV' if i == 0 else "")
		    
		    
    plot.title(f"Fake rate vs $\\eta$ for {channel}, {level}")
    plot.xlabel(r'$\eta$')
    #plot.ylabel("Factor Value")
    plot.legend()
    plot.grid(True)
    plot.savefig(f"fake_rates_lepton/plots/FR_eta_{channel}_{level}.png")

    # Plot for pt
    pt_data = level_data["pt"]
    plot.figure(figsize=(10, 6))
    #plot.plot(pt_data["values"], pt_data["factor_low_eta"], label="Low Eta", marker="o")
    #plot.plot(pt_data["values"], pt_data["factor_high_eta"], label="High Eta", marker="s")

    for i in range(len(pt_data["values"])-1):  
        plot.hlines(pt_data["factor_low_eta"][i], eta_data["values"][i], eta_data["values"][i+1], colors='blue', linewidth=2,label=r'$|\eta|<0.8$' if i == 0 else "")
    for i in range(len(pt_data["values"])-1):  
        plot.hlines(pt_data["factor_high_eta"][i], eta_data["values"][i], eta_data["values"][i+1], colors='red', linewidth=2,label=r'$|\eta|\geq0.8$' if i == 0 else "")    
    
    plot.title(f"Fake rate vs $p_{{T}}$ for {channel}, {level}")
    plot.xlabel(r'p$_{\mathrm{T}}$')
    #plot.ylabel("Factor Value")
    plot.legend()
    plot.grid(True)
    plot.savefig(f"fake_rates_lepton/plots/FR_pt_{channel}_{level}.png")

# Example usage
for ch in channels:
	for lv in levels:
		plot_factors(ch, lv)
