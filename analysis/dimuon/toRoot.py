from cafea.plotter.plotter import *
from coffea import hist, processor
import uproot3
path = 'DoubleMuon_5TeV.pkl.gz'

h = loadHistos(["DoubleMuon_dimuon.pkl.gz", "DoubleMuonLowMass_dimuon.pkl.gz"])
h = h['invmass']

out = "invmass.root"

fout = uproot3.create(out)
fout['invmass'] = hist.export1d(h)
