from config import *

#baseweb='/nfs/fanae/user/jriego/www/public/tt5TeV/'
outpath = baseweb+datatoday+'/Modeling/'
print(' >> Output = ', outpath)
if not os.path.isdir(outpath): os.makedirs(outpath)

print(' >> Loading histograms! This might take a while...')

var = 'njetsnbtags12'
plt = plotter(path, prDic=processDic, bkgList=bkglist, colors=colordic, lumi=lumi, var=var)
plt.SetLumi(lumi, "pb$^{-1}$", "5.02 TeV")
plt.SetRatio(True)
plt.SetYRatioTit('Ratio')
plt.plotData = True
plt.SetOutput(output)

# Adding hdamp and UE variations
hdampup,hdampdo = GetModSystHistos(path, '/TT_hdamp', 'hdamp', var=var) #cuidao
plt.AddExtraBkgHist([hdampup, hdampdo], add=True)
tuneup , tunedo = GetModSystHistos(path, '/TT_UE', 'UE', var=var) #cuisdao
plt.AddExtraBkgHist([tuneup, tunedo], add=True)

#'(4, 0)', '(4, 1)', '(4, $\geq$2)', '(5, 0)', '(5, 1)', '(5, $\geq$2)', '$\geq$6']
binLabels = ['(2,1)', '(2,2)', '(3,0)', '(3,1)', '(3,2)', '(4,0)', '(4,1)', '(4,2)', '(5,0)', '(5,1)', '(5,2)', '(6,*)']
 

channels = ['e', 'm']
systematics =  ['hdamp', 'UE']
level = 'incl'
for chan in channels:
    for s in systematics:
        hnom = plt.GetHistogram(var, 'tt', categories={'channel':chan, 'level':level, 'syst':'norm'})
        hup  = plt.GetHistogram(var, 'tt', categories={'channel':chan, 'level':level, 'syst':s+'Up'})
        hdo  = plt.GetHistogram(var, 'tt', categories={'channel':chan, 'level':level, 'syst':s+'Down'})
        out = outpath + 'ModUnc_ttbar_%s_%s'%(s, chan)

        fig, ax, rax = DrawComp([hnom, hup, hdo], colors=['black', 'red', 'blue'], labels=['Nominal', s+ ' Up', s+' Down'], doFill={}, outname=out, xtit='(N$_{jets}$, N$_{b-tags}$)', ytit='Events', save=False)
        # Set bin labels
        # xticks fomr 4.5 to len(binLabels)+4.5
        ax.set_xticks(np.arange(4.5, 4.5+len(binLabels), 1))
        ax.set_xticklabels(binLabels)
        rax.set_ylim(0.95, 1.05)
        # Save
        print(' >> Saving file: ', out)
        fig.savefig(out+'.pdf')
        fig.savefig(out+'.png')
