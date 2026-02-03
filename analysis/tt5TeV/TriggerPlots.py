from config import *
from cafea.plotter.plotter import GetEff

path='cafea/data/triggerSF/5TeV/triggerSFs.pkl.gz'
path='/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/triggerBiasCheck/bins_correct/'

plot = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi, var=['etatrig','etatrig_biascheck', 'pttrig', 'countstrig'])

def GetEfficiencies(var, chan, process='data'):
    hnum = plot.GetHistogram(var, categories={'level':'num', 'channel':chan, 'process':process})
    hden = plot.GetHistogram(var, categories={'level':'den', 'channel':chan, 'process':process})
    vals = GetEff(hnum, hden)
    xvals, effvals = vals
    x, xlo, xhi = xvals
    ratio, down, up = effvals
    nbins = 2
    print('hnum = ', hnum)
    print('Process: ', process, ', Channel: ', chan)
    for i in range(nbins):
        print(f'Bin {i} -- eff: {ratio[i]:.3f} +/- {up[i]:.3f}/{down[i]:.3f}')
    return ratio[:-1]
        

def DrawEffi(var, chan='m', xtit='|$\eta$|'):
    hnumMC = plot.GetHistogram(var, categories={'level':'num', 'channel':chan, 'process':'tt'})
    hdenMC = plot.GetHistogram(var, categories={'level':'den', 'channel':chan, 'process':'tt'})
    hnumDa = plot.GetHistogram(var, categories={'level':'num', 'channel':chan, 'process':'data'})
    hdenDa = plot.GetHistogram(var, categories={'level':'den', 'channel':chan, 'process':'data'})
    ynumMC = sum(hnumMC.values(overflow='all')[()])
    ydenMC = sum(hdenMC.values(overflow='all')[()])
    ynumDa = sum(hnumDa.values(overflow='all')[()])
    ydenDa = sum(hdenDa.values(overflow='all')[()])
    DrawEff(hnumMC, hdenMC, hnumDa, hdenDa, title='Efficiencies', xtit=xtit, doFill=None, mcColor='r', dataColor='k', outname=f'{var}_{chan}.png', verbose=3)
    print(f"var: {var}, chan: {chan} -- MC: {ynumMC:.1f}/{ydenMC:.1f} = {ynumMC/ydenMC:.3f}, Data: {ynumDa:.1f}/{ydenDa:.1f} = {ynumDa/ydenDa:.3f}")



def ReportTriggerEfficiencies(chan='e'):
    var = 'countstrig'
    hnumMC = plot.GetHistogram(var, categories={'level':'t', 'channel':chan, 'process':'tt'})
    hdenMC = plot.GetHistogram(var, categories={'level':'tem', 'channel':chan, 'process':'tt'})
    hnumDa = plot.GetHistogram(var, categories={'level':'num', 'channel':chan, 'process':'data'})
    hdenDa = plot.GetHistogram(var, categories={'level':'den', 'channel':chan, 'process':'data'})
    valsMC = GetEff(hnumMC, hdenMC)
    xvals, effvals = valsMC
    eff, do, up = effvals
    effMC = eff[0]; doMC = do[0]; upMC = up[0]
    valsDa = GetEff(hnumDa, hdenDa)
    xvals, effvals = valsDa
    eff, do, up = effvals
    effDa = eff[0]; doDa = do[0]; upDa = up[0]
    print('>>> Channel: ', chan, ' <<<')
    print('Efficiencies for data: %.3f +/- %.3f/%.3f' % (effDa, upDa, doDa))
    print('Efficiencies for MC  : %.3f +/- %.3f/%.3f' % (effMC, upMC, doMC))
    SF = effDa/effMC
    SFup = (effDa+upDa)/(effMC-doMC)
    SFdo = (effDa-doDa)/(effMC+upMC)
    Xmc, Ymc = GetEff(hnumMC, hdenMC)
    Xda, Yda = GetEff(hnumDa, hdenDa)
    ratio, do, up = GetRatioAssymetricUncertainties(Yda[0], Yda[1], Yda[2], Ymc[0], Ymc[1], Ymc[2])

    print('SF: %.3f +/- %.3f/%.3f' % (SF, SFup-SF, SF-SFdo))
    print('SFandrea: ',ratio, ratio+up, ratio+do)
    

    
#ReportTriggerEfficiencies('m')
#ReportTriggerEfficiencies('e')
#exit()
GetEfficiencies('etatrig', 'm', 'data')
GetEfficiencies('etatrig', 'm', 'tt')
GetEfficiencies('etatrig', 'e', 'data')
GetEfficiencies('etatrig', 'e', 'tt')

DrawEffi('etatrig', 'm')
DrawEffi('pttrig', 'm', '$p_\mathrm{T}$ (GeV)')
DrawEffi('etatrig', 'e')
DrawEffi('pttrig', 'e', '$p_\mathrm{T}$ (GeV)')

DrawEffi('etatrig_biascheck', 'm')
DrawEffi('etatrig_biascheck', 'e')
