'''
 Usage: python ljets/plotSyst.py -p histos/2jun2022_btagM_jet25/TT.pkl.gz
'''

from config import *

variation='mtop'
ch='em';  level='dilep'
plt = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi)

outpath = '/nfs/fanae/user/andreatf/www/private/ttrun3/withLepSF_withoutJECPU/syst/' 
outpath = '/nfs/fanae/user/andreatf/www/private/ttrun3/withLepSF_withoutJECPU_metfiltersOverlap_correctLepSF_recoMuonSF_PU_triggerSF/syst/' 

plt.SetOutpath(outpath)
plt.SetLumi(lumi, "pb$^{-1}$", "13.6 TeV")
plt.SetYRatioTit('Ratio')
output = "%s_%s_%s_%s"%(var, ch, level,variation)
plt.SetOutput(output)
plt.SetRatioRange(0.95,1.05)
#plt.SetRange(yRange=[0,2.5])

hdampup,hdampdo = GetModSystHistos(path, 'TTTo2L2Nu_hdamp', 'hdamp', var=var)
#tuneup , tunedo = GetModSystHistos(path, 'TT_UE', 'UE', var=var)
mtopup,mtopdo = GetModSystHistos(path, 'TTTo2L2Nu_mtop', 'mtop', var=var)

plt.AddExtraBkgHist([hdampup, hdampdo, mtopup, mtopdo], add=True)
'''
categories = {'channel':'em', 'level': 'g2jets', 'sign':'OS'}
for c in categories:
   hdampup = hdampup.integrate(c, categories[c])
   hdampdo = hdampdo.integrate(c, categories[c])
hdampup = hdampup.integrate('syst', 'hdampUp')
hdampdo = hdampdo.integrate('syst', 'hdampDown')
h = plt.GetHistogram('counts', 'tt', {'channel':'em', 'level': 'g2jets', 'sign':'OS','syst':'norm'})
#print(hdampup.values(),hdampdo.values(),h.values())
for c in categories:
   mtopup = mtopup.integrate(c, categories[c])
   mtopdo = mtopdo.integrate(c, categories[c])
mtopup = mtopup.integrate('syst', 'mtopUp')
mtopdo = mtopdo.integrate('syst', 'mtopDown')
h = plt.GetHistogram('counts', 'tt', {'channel':'em', 'level': 'g2jets', 'sign':'OS','syst':'norm'})
print(mtopup.values(),mtopdo.values(),h.values())
'''

def DrawComp(var, process, categories, labels=[], colors=[], lineStyle=[], scale=[]):
  plt.DrawComparison(var, process, categories, labels, colors, lineStyle, scale)

def DrawSyst(var, process, syst): 
  colors = ['k', 'r', 'b']
  normdict = {'channel':ch, 'level':level, 'syst':'norm'}
  h = plt.GetHistogram(var, process, normdict)
  systlist = [x.name for x in h.identifiers('syst')]
  upsyst = syst+'Up'; dosyst = syst+'Down';
  #if not upsyst in
  #dic = [{'channel':ch, 'level':level, 'syst':'norm'}, {'channel':ch, 'level':level, 'syst':syst+'Up'}, {'channel':ch, 'level':level, 'syst':syst'Do} ]

cat = [{'channel':ch, 'level':level, 'sign':'OS','syst':'norm'}, {'channel':ch, 'level':level, 'sign':'OS','syst':variation+'Up'}, {'channel':ch, 'level':level, 'sign':'OS', 'syst':variation+'Down'}]
if var is None: var = 'njets'
pr = 'tt'
labels = ['Nominal', 'up', 'down']
#labels = ['nvtx', 'calo', 'chargedPU']
colors = ['k', 'r', 'b']

DrawComp(var, pr, cat, labels, colors)

