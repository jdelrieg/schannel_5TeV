'''
 Estimate the QCD background
 Usagel: 
   qcd = QCD(path, prDic=processDic, bkglist=bkglist, lumi=lumi)
   hqcd = qcd.GetQCD(var, categories)
'''

from analysis.tt5TeV.config import *
import time

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
  if not CheckHistoCategories(plt.GetHistogram(var), {'channel' : 'e_fake', 'process':['data']+bkglist_noQCD}, checkValues=True) and not CheckHistoCategories(plt.GetHistogram(var), {'channel' : 'm_fake', 'process':['data']+bkglist_noQCD}, checkValues=True):
    print(', ', var, end='')
    continue
  variables.append(var)

# Channels and levels
print('Getting levels and systematics...')
channels    = ['e', 'm']
levels      = [x.name for x in list(plt.GetHistogram(variables[0]).identifiers('level'))]
systematics = [x.name for x in list(plt.GetHistogram(variables[0]).identifiers('syst'))]
if 'norm' in systematics: systematics.remove('norm')

print('--------------------------------------------------')
print('Variables:', variables)
print('Channels:', channels)
print('Levels:', levels)
print('Systematics:', systematics)
print('--------------------------------------------------')

def GetQCDforVar(var):
  ''' Compute (data - MC) for a given variable (dense_axis) -- all channels and levels!! '''
  catfake = {'channel' : ['e_fake', 'm_fake']}
  h_data_fake = plt.GetHistogram(var, ['data']     , catfake, keepCats=True)#.group('channel', hist.Cat("channel", "channel"), {'e_fake':'e', 'm_fake':'m'}).group('process', hist.Cat("process", "process"), {'QCD':'data'} )
  h_mc_fake   = plt.GetHistogram(var, bkglist_noQCD, catfake, keepCats=True)#.group('channel', hist.Cat("channel", "channel"), {'e_fake':'e', 'm_fake':'m'}).group('process', hist.Cat("process", "process"), {'QCD':bkglist_noQCD})
  h_data_fake = GroupKeepOrder(h_data_fake, [['channel', 'channel', {'e':'e_fake', 'm':'m_fake'}], ['process', 'process', {'QCD':'data'       }]])
  h_mc_fake   = GroupKeepOrder(h_mc_fake  , [['channel', 'channel', {'e':'e_fake', 'm':'m_fake'}], ['process', 'process', {'QCD':bkglist_noQCD}]])
  h_mc_fake.scale(-1*lumi)
  FillDataSystCategories(h_data_fake)#, var=var)
  htot = (h_data_fake + h_mc_fake)
  hd  = h_data_fake.integrate('channel', 'e').integrate('process').integrate('level', '2j1b').integrate('syst', 'norm') #4j2b
  hmc = h_mc_fake.integrate('channel', 'e').integrate('process').integrate('level', '2j1b').integrate('syst', 'norm')   #4j2b
  ht  = htot.integrate('channel', 'e').integrate('process').integrate('level', '2j1b').integrate('syst', 'norm')        #4j2b
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
  catfake = {'channel':chan+'_fake', 'level':level}
  if   sys== 1: countsData = 'counts_metl30' #l25 original #l180  ht  #l30 mtw que uso en mtw30   #l10  met con corte en 15   45 
  elif sys==-1: countsData = 'counts_metl20' #l15			#l170      #l20						  #l5						  35
  #QCD ex
  elif sys== 5: countsData = 'counts_metl29'
  elif sys== 4: countsData = 'counts_metl28'
  elif sys== 3: countsData = 'counts_metl27'
  elif sys== 2: countsData = 'counts_metl26'
  elif sys==-2: countsData = 'counts_metl24'
  elif sys==-3: countsData = 'counts_metl23'
  elif sys==-4: countsData = 'counts_metl22'
  elif sys==-5: countsData = 'counts_metl21'
    
  else:         countsData = 'counts_metl25' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_metl20    = plt.GetYields(countsData, cat    , pr='data'       , overflow='none')
  mc_metl20      = plt.GetYields(countsData, cat    , pr=bkglist_noQCD, overflow='none')
  data_metl20_fk = plt.GetYields(countsData, catfake, pr='data'       , overflow='none')
  mc_metl20_fk   = plt.GetYields(countsData, catfake, pr=bkglist_noQCD, overflow='none')
  
  if (chan=='m')&(level=='3j1b'):print(chan,level);print('data_metl20',data_metl20,'mc_metl20',mc_metl20,'data_metl20_fk',data_metl20_fk,'mc_metl20_fk',mc_metl20_fk)
  fact = (data_metl20 - mc_metl20)/(data_metl20_fk - mc_metl20_fk)
  if (chan=='m')&(level=='3j1b'):print('fact',fact,'sys',sys)
  return fact

def NormQCD(hqcd, chan, level):
  ''' Normalize QCD '''
  fact   = GetQCDnorm(chan, level)
  factUp = GetQCDnorm(chan, level, sys=1)
  factDo = GetQCDnorm(chan, level, sys=-1)
  
  factUp_5 = GetQCDnorm(chan, level, sys=5)
  factDo_5 = GetQCDnorm(chan, level, sys=-5)
  factUp_4 = GetQCDnorm(chan, level, sys=4)
  factDo_4 = GetQCDnorm(chan, level, sys=-4)
  factUp_3 = GetQCDnorm(chan, level, sys=3)
  factDo_3 = GetQCDnorm(chan, level, sys=-3)
  factUp_2 = GetQCDnorm(chan, level, sys=2)
  factDo_2 = GetQCDnorm(chan, level, sys=-2)
  if fact<0:
    fact=0.0; factUp=0.0; factDo=0.0; factUp_5=0.0; factDo_5=0.0; factUp_4=0.0; factDo_4=0.0; factUp_3=0.0; factDo_3=0.0; factUp_2=0.0; factDo_2=0.0
  cat = {'channel':chan, 'level':level}
  #for c in cat:
  #  hqcd = hqcd.integrate(c, cat[c])
  hqcd = GroupKeepOrder(hqcd, [['channel', 'channel', {cat['channel']:cat['channel']}], ['level', 'level', {cat['level']:cat['level']}]])
  hqcdUp = hqcd.copy()
  hqcdDo = hqcd.copy()
  
  hqcdUp_5 = hqcd.copy()
  hqcdDo_5 = hqcd.copy()
  hqcdUp_4 = hqcd.copy()
  hqcdDo_4 = hqcd.copy()
  hqcdUp_3 = hqcd.copy()
  hqcdDo_3 = hqcd.copy()
  hqcdUp_2 = hqcd.copy()
  hqcdDo_2 = hqcd.copy()
  #print(hqcd.sum(),hqcd.values())
  hqcd  .scale(fact)
  hqcdUp.scale(factUp)
  hqcdDo.scale(factDo)
  
  hqcdUp_5.scale(factUp_5)
  hqcdDo_5.scale(factDo_5)
  hqcdUp_4.scale(factUp_4)
  hqcdDo_4.scale(factDo_4)
  hqcdUp_3.scale(factUp_3)
  hqcdDo_3.scale(factDo_3)
  hqcdUp_2.scale(factUp_2)
  hqcdDo_2.scale(factDo_2)
  
  hqcdUp = GroupKeepOrder(hqcdUp, [['syst', 'syst', {'QCDUp':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
  hqcdDo = GroupKeepOrder(hqcdDo, [['syst', 'syst', {'QCDDown':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
  hqcd += hqcdUp
  hqcd += hqcdDo
  
  hqcdUp_5 = GroupKeepOrder(hqcdUp_5, [['syst', 'syst', {'QCDUp_5':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
  hqcdDo_5 = GroupKeepOrder(hqcdDo_5, [['syst', 'syst', {'QCDDown_5':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
  hqcd += hqcdUp_5
  hqcd += hqcdDo_5
  hqcdUp_4 = GroupKeepOrder(hqcdUp_4, [['syst', 'syst', {'QCDUp_4':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
  hqcdDo_4 = GroupKeepOrder(hqcdDo_4, [['syst', 'syst', {'QCDDown_4':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
  hqcd += hqcdUp_4
  hqcd += hqcdDo_4
  hqcdUp_3 = GroupKeepOrder(hqcdUp_3, [['syst', 'syst', {'QCDUp_3':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
  hqcdDo_3 = GroupKeepOrder(hqcdDo_3, [['syst', 'syst', {'QCDDown_3':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
  hqcd += hqcdUp_3
  hqcd += hqcdDo_3
  hqcdUp_2 = GroupKeepOrder(hqcdUp_2, [['syst', 'syst', {'QCDUp_2':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
  hqcdDo_2 = GroupKeepOrder(hqcdDo_2, [['syst', 'syst', {'QCDDown_2':'norm'}], ['process', 'process', {'QCD':'QCD'}]])
  hqcd += hqcdUp_2
  hqcd += hqcdDo_2  
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
outdict = dict(outdict)
saveHistos(path, 'QCD', outdict)

print('\nQCD saved to ', path+'/QCD.pkl.gz')

#progress100 = (progress / total) * 100
#dt = time.time() - d0
#minutes = dt / 60
#seconds = dt % 60
#expected = (dt / progress * total) if progress != 0 else 0
#expected_min = expected / 60
#expected_sec = expected % 60
#print("\r[{:<50}] {:.2f} % ({:.0f}/{:.0f}) Time: {:.0f}m {:.0f}s (expected {:.0f}m {:.0f}s)".format('#' * int(progress100/2), progress100, progress, total, minutes, seconds, expected_min, expected_sec), end="")

