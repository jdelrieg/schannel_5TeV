'''
  Usage:
    > python ttrun3/PDFscaleUncertainties.py -p histos/tttest.pkl.gz -c em -l g2jets

'''

from config import *

def Get1bPDFUnc(path, categories={}, name='PDF', sample=None, doPrint=False):
  if not 'syst' in categories: categories['syst'] = 'norm'
  if not sample is None:
    if isinstance(sample, dict):
      k = sample.keys()[0]
      categories[k] = sample[k]
    else:
      categories['sample'] = sample
  h = GetHisto(path, name, categories)
  bins, values = GetXYfromH1D(h)#, mode='edges')
  nom = values[0]
  variations = values[1:101]
  PDFunc = np.sqrt(sum( [(var - nom)*(var - nom) for var in variations] ))/nom
  alphas_up = values[102]
  alphas_do = values[101]
  AlphaSunc = abs(alphas_up - alphas_do)/2/nom
  if doPrint:
    print('PDFunc  = %1.2f %s'%(PDFunc*100, '%'))
    print('Alpha_s = %1.2f %s'%(AlphaSunc*100, '%'))
  total = np.sqrt(PDFunc*PDFunc + AlphaSunc*AlphaSunc)
  if doPrint:
    print('Total = %1.2f %s'%( total*100, ' %'))
  return total

def Get1binScaleUnc(path, categories={}, name='Scales', sample=None, doPrint=False):
  if not 'syst' in categories: categories['syst'] = 'norm'
  if not sample is None:
    if isinstance(sample, dict):
      k = sample.keys()[0]
      categories[k] = sample[k]
    else:
      categories['sample'] = sample
  h = GetHisto(path, name, categories)
  bins, values = GetXYfromH1D(h)
  muR0p5muF0p5 = values[0]
  muR0p5muF1   = values[1]
  muR0p5muF2   = values[2] # nope
  muR1muF0p5   = values[3]
  muR1muF1     = values[4] # nom
  muR1muF2     = values[5]
  muR2muF0p5   = values[6] # nope
  muR2muF1     = values[7]
  muR2muF2     = values[8]
  if doPrint:
    print('muR = 0.5, muF = 0.5 : %g   --- relative variation: %1.2f %s'%(muR0p5muF0p5, abs(muR0p5muF0p5-muR1muF1)/muR1muF1*100, '%'))
    print('muR = 0.5, muF = 1   : %g   --- relative variation: %1.2f %s'%(muR0p5muF1,   abs(muR0p5muF1-muR1muF1)/muR1muF1*100, '%' ))
    print('muR = 0.5, muF = 2   : %g   --- relative variation: %1.2f %s'%(muR0p5muF2,   abs(muR0p5muF2-muR1muF1)/muR1muF1*100, '%' ), ' --- Unphysical')
    print('muR = 1  , muF = 0.5 : %g   --- relative variation: %1.2f %s'%(muR1muF0p5,   abs(muR1muF0p5-muR1muF1)/muR1muF1*100, '%' ))
    print('muR = 1  , muF = 1   : %g   --- relative variation: %1.2f %s'%(muR1muF1,     abs(muR1muF1-muR1muF1)/muR1muF1*100, '%' ), ' --- Nominal')
    print('muR = 1  , muF = 2   : %g   --- relative variation: %1.2f %s'%(muR1muF2,     abs(muR1muF2-muR1muF1)/muR1muF1*100, '%' ))
    print('muR = 2  , muF = 0.5 : %g   --- relative variation: %1.2f %s'%(muR2muF0p5,   abs(muR2muF0p5-muR1muF1)/muR1muF1*100, '%' ), ' --- Unphysical')
    print('muR = 2  , muF = 1   : %g   --- relative variation: %1.2f %s'%(muR2muF1,     abs(muR2muF1-muR1muF1)/muR1muF1*100, '%' ))
    print('muR = 2  , muF = 2   : %g   --- relative variation: %1.2f %s'%(muR2muF2,     abs(muR2muF2-muR1muF1)/muR1muF1*100, '%' ))
  variations = [abs(muR0p5muF0p5-muR1muF1), abs(muR0p5muF1-muR1muF1), abs(muR1muF0p5-muR1muF1), abs(muR1muF2-muR1muF1), abs(muR2muF1-muR1muF1), abs(muR2muF2-muR1muF1)]
  maxvar = max(variations)
  if doPrint:
    print(' >>> Maximum variation: %1.2f %s'%(maxvar/muR1muF1*100, '%'))
  return maxvar/muR1muF1

if __name__ == '__main__':
  Get1bPDFUnc(  path, categories={'sample':'TTTo2L2Nu', 'channel':ch, 'level':level}, doPrint=True)
  Get1binScaleUnc(path, categories={'sample':'TTTo2L2Nu', 'channel':ch, 'level':level}, doPrint=True)

