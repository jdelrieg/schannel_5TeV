from config import * 

baseweb='/nfs/fanae/user/jriego/www/public/tt5TeV/'


# From master plots
fname = path + '/masterhistos/'
pmaster = plotter(fname, prDic={},  bkgList=bkglist, colors=colordic, lumi=lumi, var='master')
h = pmaster.GetHistogram('master')

#fname = path + 'Control_Plots/'
#pmaster = plotter(fname, prDic={},  bkgList=bkglist, colors=colordic, lumi=lumi, var='eeta')
#h = pmaster.GetHistogram('eeta')

dyields = {}
process = ['tt', 'tW', 'tchan', 'DY', 'WJetsH','WJetsL', 'QCD', 'data']
bkg = ['tt', 'tW', 'tchan', 'DY', 'WJetsH','WJetsL', 'QCD']
categories = [r'$e^{-}$+2j1b',r'$e^{-}$+3j1b',r'$e^{-}$+3j2b',r'$e^{+}$+2j1b',r'$e^{+}$+3j1b',r'$e^{+}$+3j2b','$\mu^{-}+2$j1b','$\mu^{-}+3$j1b','$\mu^{-}+3$j2b','$\mu^{+}+2$j1b','$\mu^{+}+3$j1b','$\mu^{+}+3$j2b']#['e+3j1b', 'e+4j1b', 'e+$\geq5$j1b', 'e+3j$\geq2$b', 'e+4j$\geq2$b', 'e+$\geq5$j$\geq2$b', '$\mu+3$j1b', '$\mu+4$j1b', '$\mu+\geq5$j1b', '$\mu+3$j$\geq2$b', '$\mu+4$j$\geq2$b', '$\mu+\geq5$j$\geq2$b']
categoriesl = ['$\ell+2$j1b','$\ell+3$j1b','$\ell+3$j2b']#'$\ell+4$j1b', '$\ell+\geq5$j1b', '$\ell+3$j$\geq2$b', '$\ell+4$j$\geq2$b', '$\ell+\geq5$j$\geq2$b']
for pr in process:
    vals = h.integrate('process', pr).integrate('syst', 'norm').values(overflow='all')[()]
    vals = vals[:-1]
    vals = vals[1:]
    vals = vals*lumi if pr != 'data' else vals
    vals = np.where(vals<0.0, 0.0, vals)
    dyields[pr] = vals

dyields['total'] = dyields['tt'] + dyields['tW'] + dyields['tchan'] + dyields['DY'] + dyields['WJetsH'] + dyields['WJetsL'] + dyields['QCD']


lyields = {}
for pr in process + ['total']:
    lyields[pr] = dyields[pr][0] + dyields[pr][1]
    
    print(pr, lyields[pr])
#datatoday="12apr23/"
output = './'
outpath = './'
print('Writing to', outpath)


def PrintYields(oname, categories, dic):
  t = OutText('./', oname, 'new', 'tex')
  ncol = len(categories)
  t.SetTexAlign('l' + ' c'*ncol)
  header = 'Yields '
  for lab in categories: header += t.vsep() + ' ' + lab
  t.line(header); t.sep()
  for pr in bkg:
    l = pr + ' '
    if pr=="WJetsH":l = "WJets (heavy)" + ' '
    if pr=="WJetsL":l = "WJets (light)" + ' '
    for v in dic[pr]:
        l += t.vsep() + ' ' + "%1.1f"%v
    t.line(l)
  t.sep()
  l = 'total '; 
  for v in dic['total']:
    l += t.vsep() + ' ' + "%1.1f"%v
  t.line(l); t.sep()
  l = 'data '
  for v in dic['data']:
    l += t.vsep() + ' ' + "%1.1f"%v
  t.line(l); t.sep()
  t.write()
  if os.path.isfile(oname+'.pdf'):
    os.system('cp %s.pdf '%oname + outpath + '/%s.pdf'%oname)


PrintYields('yields', categories, dyields)
PrintYields('yieldsl', categoriesl, lyields)


## From integrating distributions
#pp = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi)
#for chan in ['e', 'm']:
#    for lev in ['3j1b', '4j1b', 'g5j1b', '3j2b', '4j2b', 'g5j2b']:
#        pass
