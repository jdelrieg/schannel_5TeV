'''
pr['tt']['nominal']
pr['tt']['stat']
pr['tt']['norm']
pr['tt']['JES']['Up']
  Uncertainties are stored as arrays for the histograms
'''
import numpy as np
import coffea
import coffea.hist

def quad(*v):
  v = np.array(v)
  return np.sqrt( sum(v*v) )

def CheckArray(arr=None, includeStat=False):
  if isinstance(arr, dict):
    keys = list(arr.keys())
    if len(keys) == 1:
      oarr = arr[keys[0]]
      return oarr
    else:
      print("ERROR: array-like is dictionary with multiple keys: ", keys )
      return None
  elif isinstance(arr, list): 
    return np.array(arr)
  elif isinstance(arr, coffea.hist.Hist) or isinstance(arr, coffea.hist.hist_tools.Hist) or hasattr(arr, "values"):
    return CheckArray(arr.values(sumw2=includeStat))
  return arr
  
'''
def addsyst(nominal, *sys):
  # Add uncertainties for different bkg... correlated for bkg, uncorrelated for uncertainties 
  print(sys)
  diff = [CheckArray(x)-nominal for x in sys]
  vals = np.power(diff, 2)*np.sign(diff)
  up = np.sqrt(np.where(vals>0, vals, 0))
  do = np.sqrt(np.abs(np.where(vals<0, vals, 0)))
  up = np.sum(up, axis=0)
  do = np.sum(do, axis=0)
  up = nominal + up
  do = nominal - do
  print(up-nominal,nominal,nominal-do)
  return up, do
'''
def addsyst(nominal, up, do):
  ''' Add uncertainties for different bkg... correlated for bkg, uncorrelated for uncertainties '''

  diffup = [abs(CheckArray(x)-nominal) for x in up]
  diffdo = [abs(CheckArray(x)-nominal) for x in do]
  #vals = np.power(diff, 2)*np.sign(diff)
  #up = np.sqrt(np.where(vals>0, vals, 0))
  #do = np.sqrt(np.abs(np.where(vals<0, vals, 0)))

#  up = np.sum(diffup, axis=0)
  squared_diffup = np.square(diffup)
  quadratic_up = np.sum(squared_diffup, axis=0)
  
#  do = np.sum(diffdo, axis=0)
  squared_diffdo = np.square(diffdo)
  quadratic_do = np.sum(squared_diffdo, axis=0)
#  up = nominal + up
  up=nominal+np.sqrt(quadratic_up)
#  do = nominal - do
  do = nominal - np.sqrt(quadratic_do)
  return up, do
  
def addstat(nominal, *unc, relative=False):
  ''' Add stat uncertainties... uncorrelated and maybe relative '''
  if relative:
    totRel = np.sqrt(np.sum(np.power(nominal, 2)))
    up = nominal + totRel
    do = nominal - totRel
  else:
    diff = [CheckArray(x)-nominal for x in sys]
    vals = np.power(diff, 2)*np.sign(diff)
    up = np.sqrt(np.where(vals>0, vals, 0))
    do = np.sqrt(np.abs(np.where(vals<0, vals, 0)))
 
class UncertHisto:

  def __init__(self, verbose=1):
    self.process = []
    self.syst = []
    self.pr = {}
    self.verbose=verbose

  #################################################### Auxiliary
  def CheckArray(self=None, arr=None, includeStat=False):
    if isinstance(arr, dict):
      keys = list(arr.keys())
      if len(keys) == 1:
        oarr = arr[keys[0]]
        return oarr
      else:
        print("ERROR: array-like is dictionary with multiple keys: ", keys )
        return None
    elif isinstance(arr, list): 
      return np.array(arr)
    elif isinstance(arr, coffea.hist.Hist) or isinstance(arr, coffea.hist.hist_tools.Hist) or hasattr(arr, "values"):
      return self.CheckArray(arr.values(sumw2=includeStat))
    return arr

  def CheckFormatList(self, pr=None, getlist=True):
    ''' Return a well formated list of processes '''
    if pr is None: pr = self.process
    elif isinstance(pr, str) and ',' in pr: pr = pr.replace(' ', '').split(',')
    elif isinstance(pr, str): pr = [pr]
    #if not pr in self.process:
    #  print('WARNING: unknown process "%s"'%pr)
    #  return 0
    if not getlist:
      if len(pr) != 1: 
        print('[WARNING] More than one process but forced to return one...')
      return pr[0]
    return pr

  #################################################### Get uncertainties
  def GetNominal(self, process, nominal=None):
    ''' Get nominal '''
    if not nominal is None: return CheckArray(nominal)
    process = self.CheckFormatList(process)
    nom = self.pr[process[0]]['nominal']
    for pr in process[1:]:
      nom += self.pr[pr]['nominal']
    return nom

  def GetStatUnc(self, process, relative=False):
    ''' Return stat unc for a process '''
    return self.GetUncForSingleProcess(process, 'stat', relative=relative)

  def GetNormUnc(self, pr=None, relative=False):
    ''' Return norm unc for a process '''
    return self.GetUnc(process, 'norm', relative=relative)

  def GetUncForSingleProcess(self, process, unc, relative=False):
    ''' Get single unc for single process'''
    udic = self.pr[process]
    if not unc in udic.keys(): # Return nominal if unc does not exist
      if self.verbose >= 2: 
        print("[GetUnc] SOFT WARNING: uncertainty %s not found for process %s --- returning nominal"%(unc,process))
      up   = self.GetNominal(process)
      down = self.GetNominal(process)
      return up, down
    up   = udic[unc]['up']
    down = udic[unc]['down']
    return self.CheckArray(up), self.CheckArray(down)
 
  def GetUnc(self, process, unc, relative=False):
    ''' Get single unc for process or processes '''
    process = self.CheckFormatList(process)
    up, down = self.GetUncForSingleProcess(process[0], unc)
    for pr in process[1:]:
      u, d = self.GetUncForSingleProcess(pr, unc)
      up   += u
      down += d   
    return up, down

  def GetUnc_quad(self, process, unc, relative=False):
    ''' Get single unc for process or processes '''
    process = self.CheckFormatList(process)
    up, down = self.GetUncForSingleProcess(process[0], unc)
    nomtt1=self.GetUncForSingleProcess(process[0], 'norm')[0]
    nomtt2=self.GetUncForSingleProcess(process[0], 'norm')[1]
    nomtt=(nomtt1+nomtt2)/2
    nominal=np.zeros(len(nomtt))
    for pr in process:
      nominal=nominal+(self.GetUncForSingleProcess(pr, 'norm')[0]+self.GetUncForSingleProcess(pr, 'norm')[1])/2
    up=abs(up-nomtt)
    down=abs(down-nomtt)
    up_bis, down_bis = self.GetUncForSingleProcess(process[0], unc)
    for pr in process[1:]:
      u, d = self.GetUncForSingleProcess(pr, unc)
      u=abs(u-self.GetNominal(pr))
      d=abs(d-self.GetNominal(pr))
      up   = np.sqrt(up**2+u**2)
      down = np.sqrt(down**2+d**2)
    up=nominal+up
    down=nominal-down
    return up, down
    
  def GetStat(self, processes, relative=False):
    ''' Get total stat unc '''
    processes = self.CheckFormatList(processes)
    up , do = self.GetStatUnc(processes[0])
    up = up*up; do = do*do
    for pr in processes[1:]:
      u, d = self.GetStatUnc(pr)
      u = u*u; d = d*d
      up += u; do += d
    return np.sqrt(up), np.sqrt(do)

  def GetSystematic(self, process, unc=None, relative=False):
    ''' Get whatever uncertainties and processes '''
    if unc is None: unc = self.syst
    unc = self.CheckFormatList(unc)
    up = []; down = []
    for su in unc:
      u, d = self.GetUnc_quad(process, su, relative)
      up.append(self.CheckArray(u)); down.append(self.CheckArray(d))
    up, down = addsyst(self.CheckArray(self.GetNominal(process)), up, down)
    print('process',process,'\n unc',unc,'\n up',up,'\n down',down)
    return up, down

  def GetTotal(self, process, unc=None, relative=False):
    syst_up, syst_down = self.GetSystematic(process, unc, relative)
    stat_up, stat_down = self.GetStat(process, relative)
    nom = self.CheckArray(self.GetNominal(process))
    up = nom + np.sqrt( np.power(syst_up  -nom, 2) + np.power(stat_up  -nom, 2) )
    do = nom - np.sqrt( np.power(syst_down-nom, 2) + np.power(stat_down-nom, 2) )
    return up, do




  def GetDiff(self, up, down=None, nominal=None):
    ''' '''
    pass
  
 
  ########################################################## Fill and compute
  def AddProcess(self, nominal, process, stat=None, norm=None):
    self.process.append(process)
    if stat == None and hasattr(nominal, 'values'):
      nominal, stat = self.CheckArray(nominal, True)
    else: 
      nominal = self.CheckArray(nominal)
    self.pr[process] = {}
    self.pr[process]['nominal'] = nominal
    if norm is not None: self.SetNormUncForProcess(process, norm)
    if stat is not None: self.SetStatUncForProcess(process, stat)

  def SetStatUncForProcess(self, process, statUnc, relative=True):
    ''' Set stat unc for process '''
    process = self.CheckFormatList(process, 0)
    statUnc = self.CheckArray(statUnc)
    nominal = self.CheckArray(self.GetNominal(process))
    self.pr[process]['stat'] = {'up':nominal+statUnc, 'down':nominal-statUnc, 'rel':statUnc}
    return

  def SetNormUncForProcess(self, process, fact):
    ''' Set norm unc for process '''
    process = self.CheckFormatList(process)
    if len(process) == 1: process = process[0]
    else:
      for pr in process: self.SetNormUncForProcess(pr, fact)
    nominal = self.GetNominal(process)
    normup = nominal*(1+fact)
    normdo = nominal*(1-fact)
    self.pr[process]['norm'] = {'up':normup, 'down':normdo}
    return

  def SetNormUncForProcess_2(self, process, fact):
    ''' Set norm unc for process '''
    process = self.CheckFormatList(process)
    if len(process) == 1: process = process[0]
    else:
      for pr in process: self.SetNormUncForProcess(pr, fact)
    nominal = self.GetNominal(process)
    normup = nominal*(1+fact)
    normdo = nominal*(1-fact)
    self.pr[process]['norm'] = {'up':normup, 'down':normdo}
    return normup,normdo

  def AddSyst(self, name, process, up, down=None, symmetric=False, relative=False, nominal=None):
    ''' Add new uncertainty '''
    process = self.CheckFormatList(process)
    if len(process) == 1: process = process[0]
    else:
      for pr in process: self.AddSyst(name, pr, up, down, symmetric, relative, nominal)
    nominal = self.GetNominal(process, nominal)
    if down is None or symmetric: 
      up, down = self.GetSymmetric(nominal, up, down)
    if relative:
      up   = nominal + nominal*up
      down = nominal - nominal*down
    if name in self.pr[process].keys(): print('WARNING: uncertainty "%s" for process "%s" already exists! Overwriting...'%(name, process))
    self.pr[process][name] = {'up':up, 'down':down}
    return
    
  def GetSymmetric(self, nominal, up, down=None, relative = False):
    ''' Get symmetric up/down uncertainties '''
    up = self.CheckArray(up)
    if isinstance(nominal, str): nominal = self.GetNominal(nominal)
    if not down is None: # Get average uncertainties
      var1 = np.abs(nominal-up)
      var2 = np.abs(nominal-down)
      mean = (var1 + var2)/2
      up   = nominal + mean
      down = nominal - mean
      return up, down
    else:
      if not relative:
        difs = (up - nominal)
        down = nominal - difs
      else:
        down = up
    return up, down 

 
  def GetTotalStat(self):
    ''' Sum quadratically all stat unc '''
    pass

  def GetSums(self):
    ''' Sum all nominals for processes and variations '''
    pass

  def GetTotalSyst(self, pr=None, includeNorm=True):
    ''' Get total systematic by quadratically summing all systematics for a process (or for total bkg) '''
    pass


import uncertainties
def GetNSigDig(unc):
  return 1 if unc<10 else len("%i"%unc)

def GetStringNumUnc(num, unc):
  n = GetNSigDig(unc)
  if unc == 0.:
    nf = GetNSigDig(unc)
    s = "{:0.%if}"%nf
    return s.format(num), "0"
  f = uncertainties.ufloat(num, unc)
  s = "{:0.%iu}"%n
  s = s.format(f)
  nums, uncs = s.split('+/-')
  return nums, uncs
