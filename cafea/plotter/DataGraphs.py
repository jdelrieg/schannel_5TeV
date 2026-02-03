### Functions to draw nice plots from pandas dataframes using seaborn
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
from cafea.plotter.plotter import loadHistos

def BuildPandasDF(path, sampleName, branches, label=None, even=False, train_test_split=None, random_state=42):
  ''' open a pkl.gz sample and retrieve columns with names "branches", that must have the same length, and returns pandas dataset '''
  if isinstance(sampleName, str) and ',' in sampleName:
    sampleName = sampleName.replace(' ', '').split(',')

  if isinstance(sampleName, list):
    df = pd.concat([BuildPandasDF(path, s, branches, label) for s in sampleName], ignore_index=True)
    if train_test_split is not None:
      dftrain = df.sample(frac=train_test_split, random_state=random_state)
      dftest = df.drop(dftrain.index)
      return dftrain, dftest
    else:
      return df

  if isinstance(sampleName, dict):
    dataset_signal = pd.concat([BuildPandasDF(path, s, branches, 1) for s in sampleName['signal']], ignore_index=True)
    dataset_bkg    = pd.concat([BuildPandasDF(path, b, branches, 0) for b in sampleName['bkg']], ignore_index=True)
    if even:
        nmin = min([len(dataset_signal), len(dataset_bkg)])
        dataset_signal = dataset_signal.sample(nmin)
        dataset_bkg    = dataset_bkg   .sample(nmin)
    dataset = pd.concat([dataset_signal, dataset_bkg], ignore_index=True)
    if train_test_split is not None:
      dftrain = dataset.sample(frac=train_test_split, random_state=random_state)
      dftest = dataset.drop(dftrain.index)
      return dftrain, dftest
    else:
      return dataset

  path = os.path.join(path, sampleName+'.pkl.gz')
  h = loadHistos(path)
  df = pd.DataFrame()
  for name in branches:
    data = h[name].value
    df.insert(0, name, data)
  if label is not None:
    df.insert(0, 'label',  len(df)*[label])
  if train_test_split is not None:
    dftrain = df.sample(frac=train_test_split, random_state=random_state)
    dftest  = df.drop(dftrain.index)
    return dftrain, dftest
  else:
    return df

def DrawPairPlots(df, colnames, savefig=None, labname='label', **kwargs):
  ''' Draw pairplots from a pandas dataframe '''
  # do not modify colnames
  colnames = colnames[:]
  if len(colnames) > 5:
    print('WARNING: too many variables to draw pairplots, taking the first 5... Choose the most relevant ones!')
    colnames = colnames[:5]
  if not labname in colnames:
    colnames.append(labname)
  sns.set(style="ticks", color_codes=True)
  print('colnames',colnames)
  sns.pairplot(df[colnames], hue=labname, plot_kws=dict(marker="o", alpha=0.3, linewidth=0.), **kwargs)
  # size of the figure
  plt.gcf().set_size_inches(12, 10)
  if savefig is not None:
    if savefig.endswith('.png') or savefig.endswith('.pdf'): savefig = savefig[:-4]
    for end in ['.png', '.pdf']:
      plt.savefig(savefig+end)
    plt.clf()
dic={'2j1b_ht_atlas': "H'$_{T}$ (GeV)", '2j1b_mjj': "min m(j,j') (GeV)", '2j1b_medianDRjj': "median $\Delta$R(j,j')", '2j1b_mlb': "m(l,b) (GeV)", 
 '2j1b_mub': "m(u,b) (GeV)", '2j1b_absu0eta': "|$\eta_{\mathrm{u_{0}}}$|", '2j1b_deltaeta': "|$\eta_{\mathrm{u_{0}}}-\eta_{\mathrm{b_{0}}}$|", '2j1b_topmass': "m$_{\mathrm{top}}$ (GeV)", 
 '2j1b_mt': "m$_{T}$ (GeV)", '2j1b_sumAllPt': "H$_{T}$($\ell$,j) (GeV)", '2j1b_DptWub': "|$\Delta p_{\mathrm{T}}$(W,ub)| (GeV)", '2j1b_DphiWub': "|$\Delta \phi $(W,ub)|", '2j1b_Detalu': "|$\Delta\eta(\ell,\mathrm{u})$|"}
dic={'2j1b_ht_atlas': "H'$_{T}$ (GeV)", '2j1b_mjj': "min m(j,j') (GeV)", '2j1b_medianDRjj': "median $\Delta$R(j,j')", '2j1b_mlb': "m(l,b) (GeV)", 
 '2j1b_absu0eta': "|$\eta_{\mathrm{u_{0}}}$|", '2j1b_deltaeta': "|$\eta_{\mathrm{u_{0}}}-\eta_{\mathrm{b_{0}}}$|", '2j1b_topmass': "m$_{\mathrm{top}}$ (GeV)", 
 '2j1b_mt': "m$_{T}$ (GeV)",  '2j1b_DphiWub': "|$\Delta \phi $(W,ub)|", '2j1b_Detalu': "|$\Delta\eta(\ell,\mathrm{u})$|"}


def DrawHistos(df, colnames, savefig=None, labname='label', **kwargs):
    ''' Draw histograms from a pandas dataframe '''
    N = len(colnames)
    nrow = int(np.sqrt(N))
    ncol = int(np.ceil(N/nrow))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol*5, nrow*5))
    for i, name in enumerate(colnames):
        ax = axes[i//ncol, i%ncol]
        # only a line
        sns.kdeplot(df.loc[df[labname] == 0, name], ax=ax, label='Background', fill=True, linewidth=3, color='b')
        sns.kdeplot(df.loc[df[labname] == 1, name], ax=ax, label='Signal', fill=True, linewidth=3, color='r')
        ax.set_ylabel("Density", fontsize=18)
        ax.set_xlabel(dic[name],fontsize=18)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
        ax.set_title(r"$\bf{CMS}$ $\it{Simulation}$ $\it{Preliminary}$", fontsize=18, loc='left')# \it{Preliminary}$
        if name in ["2j1b_topmass"]: ax.set_xlim(0,350)
        if name in ["2j1b_DptWub"]: ax.set_xlim(0,1000)
        if name in ["2j1b_DptWub"]: ax.set_ylim(0,0.001)
        
    plt.tight_layout(pad=1.5)
    axes[0, 0].legend(fontsize=18)
    if savefig is not None:
        if savefig.endswith('.png') or savefig.endswith('.pdf'): savefig = savefig[:-4]
        for end in ['.png', '.pdf']:
          fig.savefig(savefig+end)
        plt.clf()

def DrawHistos_double(df_total,df_central, colnames, savefig=None, labname='label', **kwargs):
    ''' Draw histograms from a pandas dataframe '''
    N = len(colnames)
    nrow = int(np.sqrt(N))
    ncol = int(np.ceil(N/nrow))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol*5, nrow*5))
    for i, name in enumerate(colnames):
        ax = axes[i//ncol, i%ncol]
        # only a line
        sns.kdeplot(df_total.loc[df_total[labname] == 0, name], ax=ax, label='Bkg sides', fill=True, linewidth=3, color='b',alpha=0.3)
        sns.kdeplot(df_total.loc[df_total[labname] == 1, name], ax=ax, label='Sign sides', fill=True, linewidth=3, color='r',alpha=0.3)
        sns.kdeplot(df_central.loc[df_central[labname] == 0, name], ax=ax, label='Bkg Central', fill=True, linewidth=3, color='orange',alpha=0.1)
        sns.kdeplot(df_central.loc[df_central[labname] == 1, name], ax=ax, label='Sign Central', fill=True, linewidth=3, color='green',alpha=0.1)
        ax.set_ylabel("Density", fontsize=18)
        ax.set_xlabel(dic[name],fontsize=18)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
        ax.set_title(r"$\bf{CMS}$ $\it{Simulation}$ $\it{Preliminary}$", fontsize=18, loc='left')# \it{Preliminary}$
        if name in ["2j1b_topmass"]: ax.set_xlim(0,350)
        if name in ["2j1b_DptWub"]: ax.set_xlim(0,1000)
        if name in ["2j1b_DptWub"]: ax.set_ylim(0,0.001)
        
    plt.tight_layout(pad=1.5)
    axes[0, 0].legend(fontsize=12)
    if savefig is not None:
        if savefig.endswith('.png') or savefig.endswith('.pdf'): savefig = savefig[:-4]
        for end in ['.png', '.pdf']:
          fig.savefig(savefig+end)
        plt.clf()

def DrawBoxPlots(df, colnames, savefig=None, labname='label', **kwargs):
    ''' Draw boxplots from a pandas dataframe '''
    # if label is 0 and 1, change it to 'Background' and 'Signal'. Do not modify df!!
    df = df.copy()
    if df[labname].min() == 0 and df[labname].max() == 1:
        df[labname] = df[labname].replace({0:'Background', 1:'Signal'})
    N = len(colnames)
    ncol = int(np.sqrt(N))
    nrow = int(np.ceil(N/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol*4, nrow*4))
    for i, name in enumerate(colnames):
        ax = axes[i//ncol, i%ncol]
        sns.boxplot(x=labname, y=name, data=df, ax=ax, **kwargs)
    plt.tight_layout(pad=1.5)
    if savefig is not None:
      if savefig.endswith('.png') or savefig.endswith('.pdf'): savefig = savefig[:-4]
      for end in ['.png', '.pdf']:
        fig.savefig(savefig+end)
      plt.clf()

def DrawCorrelationMatrix(df, colnames, savefig=None, **kwargs):
    ''' Draw correlation matrix from a pandas dataframe '''
    sns.set(style="ticks", color_codes=True)
    sns.heatmap(df[colnames].corr(), **kwargs)
    plt.gcf().set_size_inches(20, 16)
    if savefig is not None:
      if savefig.endswith('.png') or savefig.endswith('.pdf'): savefig = savefig[:-4]
      for end in ['.png', '.pdf']:
        plt.savefig(savefig+end)
      plt.clf()



def _draw_as_table(df, pagesize):
    alternating_colors = [['white'] * len(df.columns), ['lightgray'] * len(df.columns)] * len(df)
    alternating_colors = alternating_colors[:len(df)]
    fig, ax = plt.subplots(figsize=pagesize)
    ax.axis('tight')
    ax.axis('off')
    ax.table(cellText=df.values, rowLabels=df.index, colLabels=df.columns, rowColours=['lightblue']*len(df), colColours=['lightblue']*len(df.columns), cellColours=alternating_colors, loc='center')
    return fig

def dataframe_to_pdf(df, filename, numpages=(1, 1), pagesize=(11, 8.5)):
  from matplotlib.backends.backend_pdf import PdfPages
  with PdfPages(filename) as pdf:
    nh, nv = numpages
    rows_per_page = len(df) // nh
    cols_per_page = len(df.columns) // nv
    for i in range(0, nh):
        for j in range(0, nv):
            page = df.iloc[(i*rows_per_page):min((i+1)*rows_per_page, len(df)),
                           (j*cols_per_page):min((j+1)*cols_per_page, len(df.columns))]
            fig = _draw_as_table(page, pagesize)
            if nh > 1 or nv > 1:
                # Add a part/page number at bottom-center of page
                fig.text(0.5, 0.5/pagesize[0], "Part-{}x{}: Page-{}".format(i+1, j+1, i*nv + j + 1), ha='center', fontsize=8)
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()
    
def DatasetSummary(df, colnames, savefig=None, labname='label', **kwargs):
    # Group by label
    df_grouped = df.groupby(labname)
    # Describe each group
    describedf = (df_grouped[colnames].describe().transpose())
    # Add a column with difference between signal and background
    describedf['diff'] = describedf[1] - describedf[0]
    # And relative difference
    describedf['rel_diff'] = describedf['diff'] / describedf[1] * 100
    # Rename columns
    describedf.columns = ['Background (B)', 'Signal (S)', 'S - B', '(S - B)/S (%)']
    if savefig is not None:
        dataframe_to_pdf(describedf, savefig, pagesize=(11, 8.5))
    return describedf


###### Decision trees
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, auc

def DrawRanking(model, df, columns, savefig=None):
    feature_importances = pd.DataFrame(model.feature_importances_,index = df[columns].columns,columns=['importance']).sort_values('importance',ascending=False)   
    print(feature_importances)
    
    name_mapping = {'2j1b_ht_atlas': "H'$_{T}$ (GeV)", '2j1b_mjj': "min m(j,j') (GeV)", '2j1b_medianDRjj': "median $\Delta$R(j,j')", '2j1b_mlb': "m(l,b) (GeV)", 
     '2j1b_mub': "m(u,b) (GeV)", '2j1b_absu0eta': "|$\eta_{\mathrm{u_{0}}}$|", '2j1b_deltaeta': "|$\eta_{\mathrm{u_{0}}}-\eta_{\mathrm{b_{0}}}$|", '2j1b_topmass': "m$_{\mathrm{top}}$ (GeV)", 
     '2j1b_mt': "m$_{T}$ (GeV)", '2j1b_sumAllPt': "H$_{T}$($\ell$,j) (GeV)", '2j1b_DptWub': "|$\Delta p_{\mathrm{T}}$(W,ub)| (GeV)", '2j1b_DphiWub': "|$\Delta \phi $(W,ub)|", '2j1b_Detalu': "|$\Delta\eta(\ell,\mathrm{u})$|"}
    feature_importances.index = feature_importances.index.map(name_mapping)

    name_mapping = {'2j1b_ht_atlas': "H'$_{T}$ (GeV)", '2j1b_mjj': "min m(j,j') (GeV)", '2j1b_medianDRjj': "median $\Delta$R(j,j')", '2j1b_mlb': "m(l,b) (GeV)", 
     '2j1b_absu0eta': "|$\eta_{\mathrm{u_{0}}}$|", '2j1b_deltaeta': "|$\eta_{\mathrm{u_{0}}}-\eta_{\mathrm{b_{0}}}$|", '2j1b_topmass': "m$_{\mathrm{top}}$ (GeV)", 
     '2j1b_mt': "m$_{T}$ (GeV)",  '2j1b_DphiWub': "|$\Delta \phi $(W,ub)|", '2j1b_Detalu': "|$\Delta\eta(\ell,\mathrm{u})$|"}
    feature_importances.index = feature_importances.index.map(name_mapping)
    
    plt.figure(figsize=(16, 16))
    sns.barplot(x=feature_importances.importance, y=feature_importances.index, hue=feature_importances.importance,  dodge=False, palette='Blues_d')
    plt.gca().get_legend().set_title('Ranking')
    # invert order of the elements in the legend
    handles, labels = plt.gca().get_legend_handles_labels()
    # transform labels into %1.3f
    labels = [float(l) for l in labels]
    labels = ['%1.2f' % l for l in labels]
    plt.gca().legend(reversed(handles), reversed(labels), loc='lower right',fontsize=36)
    plt.tight_layout(pad=1.5)
    plt.subplots_adjust(left=0.25,top=0.9,bottom=0.1) 
    plt.xlabel('Importance',fontsize=50)
    plt.ylabel('Observable',fontsize=50)
   # plt.title('Ranking of observables',fontsize=50)
    plt.text(0., -0.5, r"$\bf{CMS}$ ", fontsize=60,fontfamily='TeX Gyre Heros', horizontalalignment='left', verticalalignment='bottom')#, transform=ax.transAxes)
    plt.text(0.047, -0.5, r"$Simulation$ ", fontsize=46,fontfamily='TeX Gyre Heros', horizontalalignment='left', verticalalignment='bottom') 
    
    plt.yticks(fontsize=25)
    plt.xticks(fontsize=30)
    #plt.legend(title='Ranking', loc='lower right', labels=['{:.3f}'.format(x) for x in feature_importances.importance])
    if savefig is not None:
      if savefig.endswith('.png') or savefig.endswith('.pdf'): savefig = savefig[:-4]
      for end in ['.png', '.pdf']:
        plt.savefig(savefig+end)
      plt.clf()
    return feature_importances.importance

def ConfusionMatrix(y_true, y_pred, savefig=None):
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(10, 10))
    sns.heatmap(cm, annot=True, fmt="d")
    plt.title('Confusion matrix')
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.tight_layout(pad=1.5)
    if savefig is not None:
      if savefig.endswith('.png') or savefig.endswith('.pdf'): savefig = savefig[:-4]
      for end in ['.png', '.pdf']:
        plt.savefig(savefig+end)
      plt.clf()
    print(classification_report(y_true, y_pred))

def DrawROC(train_true, train_pred, test_true=None, test_pred=None, savefig=None):
    fpr_train, tpr_train, _ = roc_curve(train_true, train_pred)
    roc_auc_train = auc(fpr_train, tpr_train)
    plt.figure(figsize=(10, 10))
    plt.plot(fpr_train, tpr_train, label='Train (area = %0.2f)' % roc_auc_train, linewidth=3, color='darkorange')
    if test_true is not None and test_pred is not None:
        fpr_test, tpr_test, _ = roc_curve(test_true, test_pred)
        roc_auc_test = auc(fpr_test, tpr_test)
        plt.plot(fpr_test, tpr_test, label='Test (area = %0.2f)' % roc_auc_test, linewidth=3, color='mediumpurple')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.01])
    plt.xlabel('False positive rate', fontsize=22)
    plt.ylabel('True positive rate', fontsize=22)
    plt.xticks(fontsize=21)
    plt.yticks(fontsize=21)
    #plt.text(0.001, 1.015, r"$\bf{CMS}$ Preliminary", fontsize=20, horizontalalignment='left', verticalalignment='bottom')
    plt.text(0.79, 1.02, r"(5.02 TeV)", fontsize=25)
    plt.title(r"$\bf{CMS}$ $\it{Simulation}$ $\it{Preliminary}$", fontsize=25, loc='left')# \it{Preliminary}$
    #plt.title('ROC curve', fontsize=20)
    plt.legend(loc="lower right", fontsize=22)
    plt.tight_layout(pad=1.5)
    if savefig is not None:
      if savefig.endswith('.png') or savefig.endswith('.pdf'): savefig = savefig[:-4]
      for end in ['.png', '.pdf']:
        plt.savefig(savefig+end)
      plt.clf()

def DrawSigBkgHisto(train_true, train_pred, test_true=None, test_pred=None, nbins=20, savefig=None):
    plt.figure(figsize=(10, 10))
    # Normalized histograms, filled and line
    plt.hist(train_pred[train_true==1], bins=nbins, range=(0, 1), density=True, label='Signal (train)', linewidth=3, color='r', alpha=0.4)
    plt.hist(train_pred[train_true==0], bins=nbins, range=(0, 1), density=True, label='Background (train)', linewidth=3, color='b', alpha=0.4)
    if test_true is not None and test_pred is not None:
        plt.hist(test_pred[test_true==1], bins=nbins, range=(0, 1), density=True, histtype='step', label='Signal (test)', linewidth=3, color='r')#, alpha=0.9)
        plt.hist(test_pred[test_true==0], bins=nbins, range=(0, 1), density=True, histtype='step', label='Background (test)', linewidth=3, color='b')#, alpha=0.9)
    plt.xlabel('MVA score', fontsize=22)
    plt.ylabel('Density', fontsize=22)
    plt.xlim([0.0, 1.0])
    plt.xticks(fontsize=21)
    plt.yticks(fontsize=21)
    #plt.title('Score distribution', fontsize=20)
    #plt.text(0.001, 3.58, r"$\bf{CMS}$ Preliminary", fontsize=20, horizontalalignment='left', verticalalignment='bottom')
    plt.title(r"$\bf{CMS}$ $\it{Simulation}$ $\it{Preliminary}$", fontsize=25, loc='left')# \it{Preliminary}$
    plt.text(0.79, 3.83, r"(5.02 TeV)", fontsize=25)
    plt.legend(loc="best",  fontsize=22)
    plt.tight_layout(pad=1.5)
    if savefig is not None:
      if savefig.endswith('.png') or savefig.endswith('.pdf'): savefig = savefig[:-4]
      for end in ['.png', '.pdf']:
        plt.savefig(savefig+end)
      plt.clf()



