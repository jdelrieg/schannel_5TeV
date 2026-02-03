import ROOT as r
import os
import numpy as np
from copy import deepcopy


from coffea import hist
from config import *

from cafea.plotter.plotter import saveHistos, loadHistos, DrawUncPerBin, DivideHistWithErrors
from PDFscaleUncertainties import Get1binScaleUnc, Get1bPDFUnc
import matplotlib.pyplot as matplt
import matplotlib

#import tdrstyle, CMS_lumi

entry_names = ['Data', 'tchan','tbarchan','tt', 'tW', 'WJetsL','WJetsH', 'DY', 'QCD', 'unc','unc2']   #Aqui hay que meter los procesos, datos e incertidumbres (la razon por la que hay 2 incertidumbres se explica al final)
nbins=56   #Numero de bins totales en el fit para cada canal, e y mu (8 del MVA +10*8 del absu0eta)

# Create the dictionary with nbins subentries initialized to zero
dict_pre_e={entry_name: np.zeros(nbins) for entry_name in entry_names}
dict_post_e={entry_name: np.zeros(nbins) for entry_name in entry_names}

dict_pre_m={entry_name: np.zeros(nbins) for entry_name in entry_names}
dict_post_m={entry_name: np.zeros(nbins) for entry_name in entry_names} #Definimos un diccionario para cada plot




regionkk="2j1b,3j1b,3j2b"
yearkk="e_minuse_plusm_minusm_plus"
pathkk="/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/split_charge_goodJECs_mistag_comb_btagEff/temp_cards/Unblinding/e_minuse_plusm_minusm_plus"

orderedProcesses = ["QCD","DY","WJetsH", "WJetsL", "tW", "tt", "tbarchan", "tchan"]
#orderedLegendFor2Col = ["data", "vvttv", "tw", "nonworz", "ttbar", "preunc", "dy", "postunc", "ratio"]
#orderedLegendFor1Col = ["data", "tw", "ttbar", "dy", "vvttv", "nonworz", "preunc", "postunc", "ratio"]
orderedLegendFor2Col = ["data", "QCD","DY","WJetsH", "WJetsL", "tW", "tt", "tbarchan", "tchan","postunc"]
orderedLegendFor1Col = ["data", "QCD","DY","WJetsH", "WJetsL", "tW", "tt", "tbarchan", "tchan", "postunc"]

dictBinEdgesRegions = {
    "ch1"      : [0.1,0.9],
    "ch2"      : [0,5],
    "ch3"      : [0,5],
    "ch4"      : [0.1,0.9],
    "ch5"      : [0,5],
    "ch6"      : [0,5],
    "ch7"      : [0.1,0.9],
    "ch8"      : [0,5],
    "ch9"      : [0,5],
    "ch10"      : [0.1,0.9],
    "ch11"      : [0,5],
    "ch12"      : [0,5]
}

dictBinsCenterRegions = {
    "ch1"      : [0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85],
    "ch2"      : [0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75],
    "ch3"      : [0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75],
    "ch4"      : [0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85],
    "ch5"      : [0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75],
    "ch6"      : [0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75],
    "ch7"      : [0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85],
    "ch8"      : [0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75],
    "ch9"      : [0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75],
    "ch10"      : [0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85],
    "ch11"      : [0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75],
    "ch12"      : [0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75]
}




def producePlots(year, region, path):  #Esta primera parte es la herencia de CMSSW que se usa para obtener los numeros

  preFitHists = {}
  preFitHistsUnc = {}
  for key in ["prefit", "fit_s"]:
    accbin=0
    cambioch=10
    filename = path + "/fitDiagnostics{y}_{r}.root".format(y = year.replace(",", ""), r = region.replace(",", ""))
    outpath = path + "/plots_result_{y}_{r}/".format(y = year.replace(",", ""), r = region.replace(",", ""))
    #Directories with the pre or postfit results
    maindir = "shapes_%s" %(key)
    
    f = r.TFile.Open(filename)
    # Extract the fit result
    rooFit = f.Get("fit_s")
    rooFitList = rooFit.floatParsFinal()
    for param in rooFitList:
      if param.GetName() == "r":
        r_tW = param.getVal()
        break
              
    # Change to the directory containing the prefit/postfit plots
    shapes = f.Get(maindir)
    shapes.cd()
    				
	  # This directory contains a list with the different channels
    dirs = shapes.GetListOfKeys()
   # print(dirs)
    dirnames = []
    for dire in dirs:
      dirnames.append(dire.GetName())
   # print(dirnames)
    dirnames=[ 'ch1', 'ch2','ch3','ch4','ch5','ch6','ch7','ch8','ch9','ch10','ch11','ch12']  #Reordenamos los canales para que salgan como tenemos en mente
    
    for dire in dirnames:
      orderedProcesses = ["QCD","DY","WJetsH", "WJetsL", "tW", "tt", "tbarchan", "tchan"]
      #if dire=='ch6':	orderedProcesses = ["QCD", "DY", "tchan", "tW", "tt"]  
      if dire in ['ch3','ch5','ch6','ch9']:	orderedProcesses = ["QCD","WJetsH", "WJetsL", "tW", "tt", "tbarchan", "tchan"]																			 #Esto es especifico de tchannel pq con la seleccion explotaban las incertidumbres al anyadir DY en esos canales
																			     #(baja contribucion y altas uncs) En general no habria que eliminar nada aqui y con tener una unica lista de orderedprocesses deberia de valer
      if dire in ['ch1','ch4','ch7','ch10']: cambioch=8			
      if dire not in ['ch1','ch4','ch7','ch10']: cambioch=10   #cambioch es un acumulador de bines (8 para ch1,ch5,ch7,ch10 pq es el numero de bines en MVAscore distribution)
																					#y 10 para el resto pq es el numero de bines en absu0eta
      
      #TotalHisto MC
      subdir = shapes.Get(dire)
      subdir.cd()      
      htotal = subdir.Get("total")

      # Get the data points
      gr = subdir.Get("data")
      dataNpoints = gr.GetN()
      datapointsY = gr.GetY()
      datapointsX = gr.GetX()
      uncpointsHigh = gr.GetEYhigh()
      uncpointsLow = gr.GetEYlow()
      print('dataNpoints',dataNpoints,'datapointsY',np.array(datapointsY),'dir',dire)
      #if dire in ['ch1','ch4','ch7','ch10']:dataNpoints=8;datapointsY=np.array(datapointsY)[:8]
      if (key=='prefit') & (dire in ['ch1','ch2','ch3','ch4','ch5','ch6']): dict_pre_e['Data'][accbin:accbin+dataNpoints]=np.array(datapointsY)
      if (key=='prefit') & (dire in ['ch7','ch8','ch9','ch10','ch11','ch12']): dict_pre_m['Data'][accbin-nbins:accbin+dataNpoints-nbins]=np.array(datapointsY)
      if (key=='fit_s')& (dire in ['ch1','ch2','ch3','ch4','ch5','ch6']):dict_post_e['Data'][accbin:accbin+dataNpoints]=np.array(datapointsY)
      if (key=='fit_s')& (dire in ['ch7','ch8','ch9','ch10','ch11','ch12']):dict_post_m['Data'][accbin-nbins:accbin+dataNpoints-nbins]=np.array(datapointsY)    #Nos cargamos en el diccionario en la key de Data los datos para cada bin en cada canal
																																								#Observar que hacemos una por plot/diccionario
      
      subdirNew = {}
      for histoName in orderedProcesses:
        h = subdir.Get(histoName)
        htoFill = r.TH1F(histoName+"V2","", len(dictBinsCenterRegions[dire]),dictBinEdgesRegions[dire][0],dictBinEdgesRegions[dire][1])
        for bin in range(1, len(dictBinsCenterRegions[dire]) + 1):
          htoFill.SetBinContent(bin, h.GetBinContent(bin))
          #print('histoName',histoName,'getcontent',h.GetBinContent(bin),'index',accbin+bin-1)
          
          if (key=='prefit') & (dire in ['ch1','ch2','ch3','ch4','ch5','ch6']):dict_pre_e[histoName][accbin+bin-1]=h.GetBinContent(bin);    dict_pre_e['unc'][accbin+bin-1]+=h.GetBinError(bin)
          if (key=='prefit') & (dire in ['ch7','ch8','ch9','ch10','ch11','ch12']): dict_pre_m[histoName][accbin+bin-1-nbins]=h.GetBinContent(bin);  dict_pre_m['unc'][accbin+bin-1-nbins]+=h.GetBinError(bin)  #truco -nbins pq cada dict tiene nbins bines (contenedores)en total, y al tener 82 bines, si queremos separar en e y mu hay que hacer este apanyo
          if (key=='fit_s')& (dire in ['ch1','ch2','ch3','ch4','ch5','ch6']):dict_post_e[histoName][accbin+bin-1]=h.GetBinContent(bin);   dict_post_e['unc'][accbin+bin-1]+=h.GetBinError(bin)
          if (key=='fit_s')& (dire in ['ch7','ch8','ch9','ch10','ch11','ch12']):dict_post_m[histoName][accbin+bin-1-nbins]=h.GetBinContent(bin); dict_post_m['unc'][accbin+bin-1-nbins]+=h.GetBinError(bin) #Nos cargamos en el diccionario en la key de cada proceso los yields de cada proceso para cada bin en cada canal
          htoFill.SetBinError(bin, h.GetBinError(bin))																																						#Tb aqui acumulamos la unc de tipo1
          #print(dire,key,bin,histoName)
          #print(h.GetBinError(bin))
                       

      # now draw error bands
      htoFill = r.TH1F(htotal.GetName()+"V2","", len(dictBinsCenterRegions[dire]),dictBinEdgesRegions[dire][0],dictBinEdgesRegions[dire][1])
            
      for bin in range(1, len(dictBinsCenterRegions[dire]) + 1):
        htoFill.SetBinContent(bin, htotal.GetBinContent(bin))
        #print('getcontent_eror',htotal.GetBinContent(bin))
        htoFill.SetBinError(bin, htotal.GetBinError(bin))
        if (key=='prefit') & (dire in ['ch1','ch2','ch3','ch4','ch5','ch6']):    dict_pre_e['unc2'][accbin+bin-1]=htotal.GetBinError(bin)
        if (key=='prefit') & (dire in ['ch7','ch8','ch9','ch10','ch11','ch12']):  dict_pre_m['unc2'][accbin+bin-1-nbins]=htotal.GetBinError(bin)  #truco -nbins pq cada dict tiene nbins bines (contenedores)en total
        if (key=='fit_s')& (dire in ['ch1','ch2','ch3','ch4','ch5','ch6']): dict_post_e['unc2'][accbin+bin-1]=htotal.GetBinError(bin)
        if (key=='fit_s')& (dire in ['ch7','ch8','ch9','ch10','ch11','ch12']):dict_post_m['unc2'][accbin+bin-1-nbins]=htotal.GetBinError(bin)    #Aqui acumulamos la unc de tipo 2
      
      accbin+=cambioch  #Actualizamos el contador de bines sumandole cambioch (que son los bines que tiene la distribucion que toque)


  return


producePlots(yearkk, regionkk.replace(",", ""),pathkk)

############################################################################ Up to here we were using the combine tool for plots, from here to the bottom we will try to use the structure of coffea for the combined plot ################################
#		A partir de aqui se usan los numeros que tenemos almacenados en los 4 diccionarios para pintar en el estilo de coffea


hdists_pre_e   = hist.Hist("Events", hist.Cat("process", "process"), hist.Bin("shapes", "Category", 2*8 + 10*4, -0.5, (2*8 + 10*4)-0.5))     #en todo este tramo ya usamos la arquitectura de coffea (definimos coffea.hist)
hdists_pre_m   = hist.Hist("Events", hist.Cat("process", "process"), hist.Bin("shapes", "Category", 2*8 + 10*4, -0.5, (2*8 + 10*4)-0.5))     #La cosa es definir un histograma por cada uno de los 6 plots que vamos a sacar
hdists_post_e   = hist.Hist("Events", hist.Cat("process", "process"), hist.Bin("shapes", "Category", 2*8 + 10*4, -0.5, (2*8 + 10*4)-0.5))    #Ademas, por como se usa el plotter, hay que definirlos 3 veces (una para los procesos,
hdists_post_m   = hist.Hist("Events", hist.Cat("process", "process"), hist.Bin("shapes", "Category", 2*8 + 10*4, -0.5, (2*8 + 10*4)-0.5))    #otra para las incertidumbres y total [suma de los procesos] y otra para los datos)

hdistsu_pre_e   = hist.Hist("Events", hist.Cat("process", "process"), hist.Bin("shapes", "Category", 2*8 + 10*4, -0.5, (2*8 + 10*4)-0.5))    #Fijarse que el numero de bins es 58 (8 del MVA +10*5 de los medianDRjj)
hdistsu_pre_m   = hist.Hist("Events", hist.Cat("process", "process"), hist.Bin("shapes", "Category", 2*8 + 10*4, -0.5, (2*8 + 10*4)-0.5)) 
hdistsu_post_e   = hist.Hist("Events", hist.Cat("process", "process"), hist.Bin("shapes", "Category", 2*8 + 10*4, -0.5, (2*8 + 10*4)-0.5)) 
hdistsu_post_m   = hist.Hist("Events", hist.Cat("process", "process"), hist.Bin("shapes", "Category", 2*8 + 10*4, -0.5, (2*8 + 10*4)-0.5))
 
hdists1_pre_e   = hist.Hist("Events", hist.Cat("process", "process"), hist.Bin("shapes", "Category", 2*8 + 10*4, -0.5, (2*8 + 10*4)-0.5)) 
hdists1_pre_m   = hist.Hist("Events", hist.Cat("process", "process"), hist.Bin("shapes", "Category", 2*8 + 10*4, -0.5, (2*8 + 10*4)-0.5)) 
hdists1_post_e   = hist.Hist("Events", hist.Cat("process", "process"), hist.Bin("shapes", "Category", 2*8 + 10*4, -0.5, (2*8 + 10*4)-0.5)) 
hdists1_post_m   = hist.Hist("Events", hist.Cat("process", "process"), hist.Bin("shapes", "Category", 2*8 + 10*4, -0.5, (2*8 + 10*4)-0.5))   
																																	   

																																	   
random_data = np.random.normal(loc=5, scale=2, size=1000)  #Iniciailizacion de numeros random para rellenar de primeras los histogramas

for pr in orderedProcesses:hdists_pre_e.fill(shapes=random_data,process=pr);hdists_pre_m.fill(shapes=random_data,process=pr);hdists_post_e.fill(shapes=random_data,process=pr);hdists_post_m.fill(shapes=random_data,process=pr); hdists_pre_e.values()[(pr,)][:]=dict_pre_e[pr][:nbins]; hdists_pre_m.values()[(pr,)][:]=dict_pre_m[pr][:nbins]; hdists_post_e.values()[(pr,)][:]=dict_post_e[pr][:nbins]; hdists_post_m.values()[(pr,)][:]=dict_post_m[pr][:nbins]
# en la linea de arriba primero inicializo los histograms con esos numeros random y luego relleno con los valores de los diccionarios que antes fui rellenado en el produce plots


hdistsu_pre_e.fill(shapes=random_data,process='total');hdistsu_pre_m.fill(shapes=random_data,process='total');hdistsu_post_e.fill(shapes=random_data,process='total');hdistsu_post_m.fill(shapes=random_data,process='total');
hdistsu_pre_e.fill(shapes=random_data,process='unc');hdistsu_pre_m.fill(shapes=random_data,process='unc');hdistsu_post_e.fill(shapes=random_data,process='unc');hdistsu_post_m.fill(shapes=random_data,process='unc');
hdistsu_pre_e.values()[('unc',)][:]=dict_pre_e['unc'][:nbins];hdistsu_pre_m.values()[('unc',)][:]=dict_pre_m['unc'][:nbins];hdistsu_post_e.values()[('unc',)][:]=dict_post_e['unc'][:nbins];hdistsu_post_m.values()[('unc',)][:]=dict_post_m['unc'][:nbins];

hdistsu_pre_e.fill(shapes=random_data,process='unc2');hdistsu_pre_m.fill(shapes=random_data,process='unc2');hdistsu_post_e.fill(shapes=random_data,process='unc2');hdistsu_post_m.fill(shapes=random_data,process='unc2');
hdistsu_pre_e.values()[('unc2',)][:]=dict_pre_e['unc2'][:nbins];hdistsu_pre_m.values()[('unc2',)][:]=dict_pre_m['unc2'][:nbins];hdistsu_post_e.values()[('unc2',)][:]=dict_post_e['unc2'][:nbins];hdistsu_post_m.values()[('unc2',)][:]=dict_post_m['unc2'][:nbins];
#Hago lo mismo que en los hists de procesos para las incertidumbres 


hdistsu_pre_e.values()[('total',)][:]=dict_pre_e['tchan'][:nbins]+dict_pre_e['tbarchan'][:nbins]+dict_pre_e['tt'][:nbins]+dict_pre_e['tW'][:nbins]+dict_pre_e['WJetsL'][:nbins]+dict_pre_e['WJetsH'][:nbins]+dict_pre_e['QCD'][:nbins]+dict_pre_e['DY'][:nbins]
hdistsu_pre_m.values()[('total',)][:]=dict_pre_m['tchan'][:nbins]+dict_pre_m['tbarchan'][:nbins]+dict_pre_m['tt'][:nbins]+dict_pre_m['tW'][:nbins]+dict_pre_m['WJetsL'][:nbins]+dict_pre_m['WJetsH'][:nbins]+dict_pre_m['QCD'][:nbins]+dict_pre_m['DY'][:nbins]
hdistsu_post_e.values()[('total',)][:]=dict_post_e['tchan'][:nbins]+dict_post_e['tbarchan'][:nbins]+dict_post_e['tt'][:nbins]+dict_post_e['tW'][:nbins]+dict_post_e['WJetsL'][:nbins]+dict_post_e['WJetsH'][:nbins]+dict_post_e['QCD'][:nbins]+dict_post_e['DY'][:nbins]
hdistsu_post_m.values()[('total',)][:]=dict_post_m['tchan'][:nbins]+dict_post_m['tbarchan'][:nbins]+dict_post_m['tt'][:nbins]+dict_post_m['tW'][:nbins]+dict_post_m['WJetsL'][:nbins]+dict_post_m['WJetsH'][:nbins]+dict_post_m['QCD'][:nbins]+dict_post_m['DY'][:nbins] #Relleno el total sumando las contribuciones de cada proceso

hdists1_pre_e.fill(shapes=random_data,process='Data')
hdists1_pre_e.values()[('Data',)][:]=dict_pre_e['Data'][:nbins]
hdists1_pre_m.fill(shapes=random_data,process='Data')
hdists1_pre_m.values()[('Data',)][:]=dict_pre_m['Data'][:nbins]
hdists1_post_e.fill(shapes=random_data,process='Data')
hdists1_post_e.values()[('Data',)][:]=dict_post_e['Data'][:nbins]
hdists1_post_m.fill(shapes=random_data,process='Data')
hdists1_post_m.values()[('Data',)][:]=dict_post_m['Data'][:nbins]   #Misma historia con datos


print(dict_post_e)

def DrawShapesMaser(fname, chan='m', doData=True,channel='epre',chooseunc='type1'):  #Metodo para hacer el plot (fundamentalmente lo que se hace esta en el Stack_2 del plotter.py)
    ''' Drwa the master histogram '''												#Aspectos importantes: channel epre,epost,mpre y mpost para los 4 plots; chooseunc: tipo 1 o tipo 2
    outpath= pathkk																	#Y en Stack_2 en funcion del channel se carga 1 de los 4 sets de histogramas (normal, u y data)
    if not os.path.exists(outpath):													#El resto es sobretodo estetica.
        os.makedirs(outpath)
    outname = 'shapes_'+channel
    
    plt = plotter(fname, prDic={},  bkgList=['tchan', 'tbarchan','tt', 'tW','WJetsL', 'WJetsH',  'DY', 'QCD'], colors=colordic, lumi=lumi, var='shapes')
    plt.SetLumi(lumi, "pb$^{-1}$", "5.02 TeV")
    plt.SetCategories({"channel":chan})
    plt.plotData = doData
    if channel=="epre":hdists=hdists_pre_e;hdists1=hdists1_pre_e;hdistsu=hdistsu_pre_e;
    if channel=="mpre":hdists=hdists_pre_m;hdists1=hdists1_pre_m;hdistsu=hdistsu_pre_m;
    if channel=="epost":hdists=hdists_post_e;hdists1=hdists1_post_e;hdistsu=hdistsu_post_e;
    if channel=="mpost":hdists=hdists_post_m;hdists1=hdists1_post_m;hdistsu=hdistsu_post_m;
    fig, ax, rax = plt.Stack_2(hdists,hdists1,hdistsu, xtit='', ytit='Events', dosyst=True, verbose=1, doNotSave=True,chooseunc=chooseunc)
    ax.set_ylim(0,125)        
    if chan=='e': ax.set_ylim(0,80)
    # Legend 
    handles, labels = ax.get_legend_handles_labels()

    for k, lab in diclegendlabels.items():
        if k in labels:
            labels[labels.index(k)] = lab
            
    unc_handle = mpatches.Rectangle((0, 0), 1, 1,facecolor='white', hatch="\/\/",label='Unc.',edgecolor='gray',linewidth=0)
    all_handles =  [unc_handle]+handles
    all_labels =  ['Unc.']+labels
    all_handles.pop(7)
    all_labels.pop(7)
    ax.legend(all_handles[::-1], all_labels[::-1], loc='upper right', ncol=3, fontsize=14,frameon=True,framealpha=1,facecolor='white',fancybox=False,edgecolor='white')
    if doData:
      rax.set_ylim(0.2, 1.8)
      rax.set_xticks([])
    
    sep = [7.5,17.5,27.5,35.5,45.5]; #points for vertical lines
    for axx in [ax, rax]:
      if axx is None: continue
      for x in sep:
        if axx==ax:axx.axvline(x=x, color='k', linestyle='--')#,ymax=0.65)
        if axx==rax:axx.axvline(x=x, color='k', linestyle='--')
        
    fig.subplots_adjust(right=0.97, top=0.94, bottom=0.06, left=0.08)
    fig.set_size_inches(12, 9)
    if doData==False: fig.tight_layout(rect=(0,0.02,1,1))
    ax.text(0.38, 0.93, '$e+$jets, Post-fit' if chan=='e' else '$\mu+$jets, Post-fit', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=22, color='k')
    catlabels = ['2j1b', '3j1b', '3j2b', '2j1b','3j1b','3j2b']
    xpos = [0.07, 0.22, 0.4,0.57,0.75, 0.92]
    for cat, x in zip(catlabels, xpos):
      ax.text(x, 0.65, cat, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=22, color='k')
      
    
      
    if doData:
      rax.text(0.12,  -0.26, "MVA Score", horizontalalignment='center', verticalalignment='center', transform=rax.transAxes, fontsize=20, color='k')

      ax.set_xticks(ticks=[-0.5,1.5,3.5,6.5,8,9.75,11.5,13.5,15,16.75,18.5,20.5,22,23.75,25.5],labels=[0.2,0.4,0.6,1.5,2,2.5,3,1.5,2,2.5,3,1.5,2,2.5,3],fontsize=18)
      
      rax.set_xlabel("", fontsize=20)
      rax.text(0.35,  -0.26, r"$\Delta R_\mathrm{med}$(j,j')", horizontalalignment='center', verticalalignment='center', transform=rax.transAxes, fontsize=20, color='k')
      rax.text(0.61,  -0.26, r"$\Delta R_\mathrm{med}$(j,j')", horizontalalignment='center', verticalalignment='center', transform=rax.transAxes, fontsize=20, color='k')
      rax.text(0.87,  -0.26, r"$\Delta R_\mathrm{med}$(j,j')", horizontalalignment='center', verticalalignment='center', transform=rax.transAxes, fontsize=20, color='k')
      rax.set_ylabel('Data / Fit')
      fig.set_size_inches(14, 10)
      fig.subplots_adjust(bottom=0.07) 


      for sufix in ['png', 'pdf']: fig.savefig(outpath + outname + '.' + sufix)
    else:
      ax.text(0.07,  -0.065, "MVA Score", horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20, color='k')
      ax.set_xlabel(r'$\Delta R_\mathrm{med}$(j,j)', fontsize=20)
      ax.set_xticks(ticks=[-0.5,1.5,3.5,6.5,8,9.75,11.5,13.5,15,16.75,18.5,20.5,22,23.75,25.5,27.5,29,30.75,32.5,34.5,36,37.75,39.5],labels=[0.05,0.4,0.6,1.5,2,2.5,3,1.5,2,2.5,3,1.5,2,2.5,3,1.5,2,2.5,3,1.5,2,2.5,3],fontsize=18)
      fig.set_size_inches(14, 10)
      for sufix in ['png', 'pdf']: fig.savefig(outpath + outname + '_blind.' + sufix)
    print('Saved to: ', outpath + outname + '.png')
    
fname='/nfs/fanae/user/jriego/tchannel5TeV/splitting_tchan_tbar/split_charge_goodJECs_mistag_comb_btagEff/masterhistos/master.pkl.gz'


#DrawShapesMaser(fname, 'e', True,channel='epre',chooseunc='type2')   #Y aqui es donde llamamos al metodo definido para hacer cada uno de los plots. Conviene explicar aqui lo de las uncs:
#DrawShapesMaser(fname, 'm', True,channel='mpre',chooseunc='type2')   #Type 2 es lo que nos daba el plot inicial de cmssw, que es una incertidumbre que se carga en l161 (asociada al histograma total)
DrawShapesMaser(fname, 'e', True,channel='epost',chooseunc='type2')  #Por defecto pongo esa ya que es la que presentabamos en los plots de antes.  
#DrawShapesMaser(fname, 'm', True,channel='mpost',chooseunc='type2')  #Type 1 es ir sumando la incertidumbre proceso a proceso y no tomando una total. Da lugar a una incertidumbre mayor, y en verdad creo que es mas correcto  
																	 #pero puesto que la type 2 es la que se usaba por defecto, de momento lo dejo asi
