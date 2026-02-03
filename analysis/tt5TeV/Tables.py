from config import *
import uproot
from QCD import GetQCDbkg

def GetYieldsDic(fname, listOfSyst=None, processes=None):
  ''' Get a dictionary with yields taking a combine-format rootfile and a list of systematics '''
  f = uproot.open(fname)
  if processes is None: processes = ['tt', 'tW', 'tchan', 'WJets', 'DY', 'QCD']
  yields = {}
  yields['data'] = {'': {'values': f['data_obs'].values(), 'errors': np.sqrt(f['data_obs'].values())}}
  if listOfSyst is None: listOfSyst=['']
  elif not '' in listOfSyst: listOfSyst = [''] + listOfSyst
  for pr in processes:
    yields[pr] = {}
    for syst in listOfSyst:
      for var in ['Up', 'Down']:
        if syst == '' and var == 'Down': continue
        hname = pr if syst=='' else pr+'_'+syst+var
        chsyst = syst+var if not syst=='' else ''
        yields[pr][chsyst] = {'values': f[hname].values(), 'errors':f[hname].errors()}
  return yields

binNames = ['e 0b', 'e 1b', 'e $\geq2$b', '$\mu$ 0b', '$\mu$ 1b', '$\mu$ $\geq2$b']

tdummy = OutText('', 'Yields', 'new', 'tex')
def NLine(first, labList, sep=tdummy.vsep(), form='%1.2f'):
  for l in labList:
    first += sep + ' ' + form%l
  return first

def PrintYields(path, outpath=None, outname='Yields'):
  if outpath is None: outpath = outpatho
  yields = GetYieldsDic(path, [''])
  t = OutText(outpath, outname, 'new', 'tex')

  Y = lambda pr : np.array(yields[pr]['']['values'])
  S = lambda sys : ( abs(np.array(yields[pr][sys+'Up']['values']) - Y('tt'))/Y('tt') + abs(np.array(yields[pr][sys+'Down']['values'])-Y('tt'))/Y('tt') )/2*100
  totBkg = Y('tt') + Y('tW') + Y('tchan') + Y('WJets') + Y('DY') + Y('QCD')

  # Calculate yield table
  bkgList = ['tW', 'tchan', 'WJets', 'QCD', 'DY']
  t.bar()
  ncol = 6
  t.SetTexAlign('l' + ' c'*ncol)
  header = 'Yields '
  for lab in binNames: header += t.vsep() + ' ' + lab
  t.line(header); t.sep()
  t.line( NLine('tW', Y('tW')) )
  t.line( NLine('W+jets', Y('WJets')) )
  t.line( NLine('QCD', Y('QCD')) )
  t.line( NLine('Drell--Yan', Y('DY')) )
  t.line( NLine('$\mathrm{t}\\bar{\mathrm{t}}$', Y('tt')) )
  t.sep()
  t.line( NLine('Total pred. ', totBkg) )
  t.sep()
  t.line( NLine('Data', Y('data')) )
  t.sep()
  t.write()

def CreateSystTable(path , outpath=None, systematics=['lepSF', 'btagSF', 'trigSF', 'ISR', 'FSR'], pr='tt',systNom = ''):
  if outpath is None: outpath = outpatho
  yields = GetYieldsDic(path, systematics, [pr])
  Y = lambda pr : np.array(yields[pr]['']['values'])
  S = lambda sys : ( abs(np.array(yields[pr][sys+'Up']['values']) - Y('tt'))/Y('tt') + abs(np.array(yields[pr][sys+'Down']['values'])-Y('tt'))/Y('tt') )/2*100
  # Calculate systematic table
  t = OutText(outpath, 'Systematics', 'new', 'tex')
  t.bar()
  ncol = 6
  t.SetTexAlign('l' + ' c'*ncol)
  header = 'Source '
  for lab in binNames: header += t.vsep() + ' ' + lab
  t.line(header); t.sep()
  for s in systematics:
    t.line( NLine(s, S(s)) )
  t.sep()
  t.write()

if __name__ == '__main__':
  CreateSystTable(path)
  PrintYields(path)


