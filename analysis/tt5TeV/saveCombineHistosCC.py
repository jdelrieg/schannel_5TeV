# Create a combine-like rootfile with the proper format
# Usage:
#    >> python SaveCombineHistosCC.py -p histos/pklFiles/
#    >> python SaveCombineHistosCC.py -p ../histos/20dec2021/ -o /nfs/fanae/user/juanr/CMSSW_10_2_13/src/tt5TeV/ljets --output NBtags.root


from config import *
from QCD import GetQCDbkg

def SaveCombine(verbose=False):
  plt = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi)
  outname = output if output is not None else 'temp'
  yields = {}
  out = outpatho + outname + ('.root' if not outname.endswith('.root') else '')
  if os.path.isfile(out): os.system('mv %s %s_old'%(out, out))
  fout = uproot3.create(out)
  for pr in ['tt', 'tW', 'tchan', 'WJets', 'DY', 'data', 'QCD']:
    yields[pr] = {}
    for syst in ['norm', 'lepSFUp',"lepSFDown","btagSFUp", "btagSFDown", "trigSFUp", "trigSFDown", 'ISRUp', 'ISRDown', 'FSRUp', 'FSRDown']:
      yields[pr][syst] = {}
      if syst in ['ISRUp', 'ISRDown', 'FSRUp', 'FSRDown'] and pr != 'tt': continue
      if syst != 'norm' and pr in ['data', 'QCD']: continue
      v = []; ve = []
      for ch in ['e', 'm']:
        for l in ['0b', '1b', '2b']:
          categories = {'channel':ch, 'level':l, 'syst':syst}
          if pr == 'QCD' and syst == 'norm':
            hQCD = GetQCDbkg('counts', categories)
            y = ((hQCD.values(sumw2=True))[('QCD',)])[0]
            ye= ((hQCD.values(sumw2=True))[('QCD',)])[1]
            y  = float(y ) if  y  > 0 else  0.0
            ye = float(ye) if  ye > 0 else  0.0
          elif pr == 'tt' and syst in ['hdampUp', 'hdampDown', 'UEUp', 'UEDown']:
            continue
          else:
            y, ye = plt.GetYields('counts', categories, doErr=True)
            y = y[pr]; ye = ye[pr]
            if pr == 'data': ye = np.sqrt(y)
          if verbose: print('(%s) = %1.2f +/- %1.2f'%(pr, y, ye))
          v.append(y)
          ve.append(ye)
      y = np.array(v)
      yerr = np.array(ve)
      bins = np.array(range(len(y)+1))
      yields[pr][syst]['values'] = y
      yields[pr][syst]['errors'] = yerr
      h = GetH1DfromXY( bins, y, yerr )
      prname = pr if not pr=='data' else 'data_obs'
      if pr != 'data' and syst != 'norm': prname = pr+'_'+syst
      fout[prname] = hist.export1d(h)
  print('Created file: ', out)
  fout.close()

if __name__ == '__main__':
  SaveCombine()
