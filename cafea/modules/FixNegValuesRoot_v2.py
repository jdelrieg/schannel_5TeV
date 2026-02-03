import ROOT as r
import sys, os
from copy import deepcopy


def SetNegValToZero(h):
  for i in range(h.GetNbinsX()+2):
    x = h.GetBinContent(i)
    if i < 0: h.SetBinContent(i, defaultZeroVal)
  return h


def FixNegValuesRoot(path, defaultZeroVal, verbose = False):
  if verbose: print(f"> Checking for negative or zero bins for file {path}.")
  f = r.TFile(path, 'READ')
  outhistos = []
  hnames = [k.GetName() for k in f.GetListOfKeys()]
  nFixed = 0
  for name in hnames:
    outhistos.append(deepcopy(f.Get(name).Clone()))
    fixed = 0
    for i in range(outhistos[-1].GetNbinsX()+2):
      x = outhistos[-1].GetBinContent(i)
      if x < defaultZeroVal:
        fixed = 1
        outhistos[-1].SetBinContent(i, defaultZeroVal)
    nFixed += fixed
  f.Close(); del f
  if verbose: print(f"\t- Fixed {nFixed} out of {len(hnames)} histograms.")
  os.system("mv " + path + " " + path.replace(".root", "_prefix.root"))
  of = r.TFile(path, 'RECREATE')
  for iH in outhistos: iH.Write()
  of.Close(); del of
  return


if __name__ == '__main__':
  args = sys.argv[1:]
  if len(args) == 0:
    print("Usage:\n  >> python cafea/modules/FixNegValuesRoot.py /path/to/rootfile.root")
    exit()
  fname = args[-1]
  FixNegValuesRoot(fname, verbose=True)
