from ROOT import *
gROOT.SetBatch(1)
fname = 'invmass.root'
pname = 'invmass'

outpath = '~/www/public/'

def DrawLabelName(name, X, Y, s = 0.045):
 ts = TLatex(-20.,50., name);
 ts.SetNDC();
 ts.SetTextAlign(12);
 ts.SetX(X);
 ts.SetY(Y);
 ts.SetTextFont(42);
 ts.SetTextSize(s);
 ts.SetTextSizePixels(22);
 return ts

f = TFile.Open(fname)
h = f.Get(pname)

h.SetStats(0)
h.SetTitle('')

h.GetXaxis().SetTitle('M_{#mu#mu} (GeV)')
h.GetXaxis().SetTitleSize(0.05)
h.GetXaxis().SetTitleOffset(1.12)

h.GetYaxis().SetTitle('Events / GeV')
h.GetYaxis().SetTitleSize(0.05)
h.GetYaxis().SetTitleOffset(1.12)
h.GetXaxis().SetRangeUser(0.25, 150)

nbins = h.GetNbinsX()

for i in range(nbins):
  c = h.GetBinContent(i)
  w = h.GetBinWidth(i)
  h.SetBinContent(i, c/w)

h.SetMaximum(10e8)
h.SetMinimum(1)

c = TCanvas('c', 'c', 10, 10, 800, 600)
gPad.SetTickx(); gPad.SetTicky();
c.SetLogy(); c.SetLogx()
c.SetTopMargin(0.07); c.SetBottomMargin(0.15); c.SetRightMargin(0.03); c.SetLeftMargin(0.12)

h.Draw('lc')
tex=[]
tex.append(DrawLabelName("296.1 pb^{-1} (5.02 TeV)", 0.655, 0.96, 0.048))
tex.append(DrawLabelName("#eta", 0.215, 0.74))
tex.append(DrawLabelName("#rho/#omega", 0.245, 0.77))
tex.append(DrawLabelName("#phi", 0.3, 0.773))
tex.append(DrawLabelName("J/#psi", 0.43, 0.86))
tex.append(DrawLabelName("#psi'", 0.47, 0.735))
tex.append(DrawLabelName("#Upsilon", 0.60, 0.74))
tex.append(DrawLabelName("Z", 0.893, 0.55))

for l in tex: l.Draw()

c.Print(outpath + 'InvMass_test.png', 'png')
c.Print(outpath + 'InvMass_test.pdf', 'pdf')

