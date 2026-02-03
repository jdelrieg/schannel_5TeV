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
  if not CheckHistoCategories(plt.GetHistogram(var), {'channel' : 'e_fake_plus', 'process':['data']+bkglist_noQCD}, checkValues=True) and not CheckHistoCategories(plt.GetHistogram(var), {'channel' : 'm_fake_plus', 'process':['data']+bkglist_noQCD}, checkValues=True) and not CheckHistoCategories(plt.GetHistogram(var), {'channel' : 'm_fake_minus', 'process':['data']+bkglist_noQCD}, checkValues=True) and not CheckHistoCategories(plt.GetHistogram(var), {'channel' : 'e_fake_minus', 'process':['data']+bkglist_noQCD}, checkValues=True):
    print(', ', var, end='')
    continue
  variables.append(var)

# Channels and levels
print('Getting levels and systematics...')
channels    = ['e_plus','e_minus', 'm_plus','m_minus']
levels      = [x.name for x in list(plt.GetHistogram(variables[0]).identifiers('level'))]
systematics = [x.name for x in list(plt.GetHistogram(variables[0]).identifiers('syst'))]
if 'norm' in systematics: systematics.remove('norm')
levels=['3j2b']
systematics=['norm']

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
  if (sys==0):print('fake rate',chan,level,fact);print('data_metl20',data_metl20,'mc_metl20',mc_metl20,'data_metl20_fk',data_metl20_fk,'mc_metl20_fk',mc_metl20_fk)
  #if sys==1:print('fake rate up',chan,level,fact)
  #if sys==-1:print('fake rate down',chan,level,fact)

  return fact

def GetQCDnorm_v2(chan, level, sys=0):
  cat     = {'channel':chan, 'level':level}
  catfake = {'channel':chan[0] + '_fake_' + chan[2:], 'level':level}
  countsData = 'counts' #l20          #l175      #l25 						  #l75(7.5)					  40
  data_fk = plt.GetYields(countsData, catfake, pr='data'       , overflow='none')
  mc_fk   = plt.GetYields(countsData, catfake, pr=bkglist_noQCD, overflow='none')
  
  #if sys==0:print(chan,level);print('data_fk',data_fk,'mc_fk',mc_fk)
  fact = (data_fk - mc_fk)
  if sys==0:print('fact counts',chan,level,fact)
  #if (sys==0)&(level=='1j1b'):print('fact counts',chan,level,fact);print('data_fk',data_fk,'mc_fk',mc_fk)
  return fact


def NormQCD(hqcd, chan, level):
  ''' Normalize QCD '''
  fact   = GetQCDnorm(chan, level)
  otro_factor=GetQCDnorm_v2(chan,level)
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

