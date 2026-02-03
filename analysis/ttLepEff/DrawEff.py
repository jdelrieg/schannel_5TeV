from config import *
import copy

MCprompt = ["tt", "tW", "DY", "Diboson"]
MCnonPrompt = ["WJets", "semilep"]

# Axes
# process: tt, tW, semilep, WJets, DY, Diboson
# level: dilep, g2jets, g1btag
# channel: e, m, em
# sign: OS, SS
# lepton: pass, fail, genmatch


############### Efficiencies 
def GetSSOSfact(plt, var, cat={}):
  ''' Get OS/SS factor using MC nonprompt events from WJets and tt semilep '''
  catcopy = cat.copy()
  catcopy["process"] = MCnonPrompt
  catcopy['sign'] = 'SS'
  hSS = plt.GetHistogram(var, categories=catcopy).copy()
  catcopy['sign'] = 'OS'
  hOS = plt.GetHistogram(var, categories=catcopy).copy()
  # Integrate the rest of axes
  for ax in list(hSS.axes()):
    hSS = hSS.integrate(ax)
    hOS = hOS.integrate(ax)
  fact = hOS.values()[()]/hSS.values()[()]
  return fact
  
def GetNonprompt(plt, var, cat={}):
  ''' Get nonprompt estimate using SS data '''
  catdata = cat.copy()
  catdata["sign"] = "SS"
  catMC   = cat.copy()
  catMC["sign"] = "SS"
  hdataSS = plt.GetHistogram(var, categories=catdata).group("process", hist.Cat("process", "process"), {'Nonprompt':'data'})
  hMCSS   = plt.GetHistogram(var, categories=catMC).group("process", hist.Cat("process", "process"), {'Nonprompt':MCprompt})
  rOSSS = GetSSOSfact(plt, var, cat)
  hMCSS.scale(-1*plt.lumi)
  h = (hdataSS + hMCSS)
  h.scale(rOSSS)
  return h

#################################
# Data eff = (data pass - nonprompt pass) / ( (data pass + data fail) - (nonprompt pass - nonprompt fail))
# MC   eff = (MC prompt pass) / (MC prompt pass + MC prompt fail)
# DY   eff = (MC prompt pass) / (MC prompt pass + MC prompt fail) * SF from DY

def GetDataEff(plt, var, cat={}):
  hNonprompt = GetNonprompt(plt, var, cat)
  hdata = plt.GetHistogram(var, cat)
  hdata = hdata.integrate("process", "data")
  hdataTot = hdata.sum("lepton")
  hdataPass = hdata.integrate("lepton", "pass")
  hNonpromptTot = hNonprompt.sum("lepton")
  hNonpromptPass = hNonprompt.integrate("lepton", "pass")
  hNum = hdataPass - hNonpromptPass
  hDen = hdataTot - hNonpromptTot
  eff = hNum/hDen
  return eff

def GetMCEff(plt, var, cat={}):
  hMC = plt.GetHistogram(var, cat)
  if 'procses' in list(hMC.axes()): hMC = hMC.integrate("process", MCprompt)
  hMCPass = hMC.integrate("lepton", "pass")
  hMCTot = hMCPass + hMC.integrate("lepton", "fail")
  eff = hMCPass/hMCTot
  return eff

def GetMCcorrEff(plt, var, cat={}):
  hMC = plt.GetHistogram(var, cat)
  if 'process' in list(hMC.axes()): hMC = hMC.integrate("process", MCprompt)
  hMCPass = hMC.integrate("lepton", "pass")
  hMCTot = hMCPass + hMC.integrate("lepton", "fail")
  eff = hMCPass/hMCTot
  return eff

if __name__ == "__main__":
  plt = plotter(var, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi, var=var)
