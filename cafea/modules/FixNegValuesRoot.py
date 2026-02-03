import ROOT
import uproot
import sys, os

def SetNegValToZero(h):
  for i in range(h.GetNbinsX()+2):
    x = h.GetBinContent(i)
    if i < 0: h.SetBinContent(i, 0)
  return h

def FixNegValuesRoot(path, verbose=False):
  f = ROOT.TFile.Open(path, 'update')
  hnames = [k.GetName() for k in f.GetListOfKeys()]
  nFixed = 0
  for name in hnames:
    h = f.Get(name)
    fixed = 0
    for i in range(h.GetNbinsX()+2):
      x = h.GetBinContent(i)
      if x < 0: 
        fixed = 1
        h.SetBinContent(i, 0)
    nFixed += fixed
    h.Write()
  if verbose: print(f" >> Fixed {nFixed} out of {len(hnames)} histograms.")  
  f.Close()

if __name__ == '__main__':
  args = sys.argv[1:]
  if len(args) == 0:
    print("Usage:\n  >> python cafea/modules/FixNegValuesRoot.py /path/to/rootfile.root")
    exit()
  fname = args[-1]
  FixNegValuesRoot(fname, verbose=True)
