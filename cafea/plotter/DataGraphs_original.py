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
  sns.pairplot(df[colnames], hue=labname, plot_kws=dict(marker="o", alpha=0.3, linewidth=0.), **kwargs)
  # size of the figure
  plt.gcf().set_size_inches(12, 10)
  if savefig is not None:
    if savefig.endswith('.png') or savefig.endswith('.pdf'): savefig = savefig[:-4]
    for end in ['.png', '.pdf']:
      plt.savefig(savefig+end)
    plt.clf()

def DrawHistos(df, colnames, savefig=None, labname='label', **kwargs):
    ''' Draw histograms from a pandas dataframe '''
    N = len(colnames)
    nrow = int(np.sqrt(N))
    ncol = int(np.ceil(N/nrow))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol*4, nrow*4))
    for i, name in enumerate(colnames):
        ax = axes[i//ncol, i%ncol]
        # only a line
        sns.kdeplot(df.loc[df[labname] == 0, name], ax=ax, label='Background', fill=True, linewidth=3, color='b')
        sns.kdeplot(df.loc[df[labname] == 1, name], ax=ax, label='Signal', fill=True, linewidth=3, color='darkorange')
    plt.tight_layout(pad=1.5)
    axes[0, 0].legend()
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
    custom_palette = sns.color_palette("YlOrBr", as_cmap=True)
    sns.heatmap(df[colnames].corr(),cmap=custom_palette ,**kwargs)
    plt.gcf().set_size_inches(12, 10)
    plt.tight_layout(rect=(0,0.02,1,1))
 
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
    plt.figure(figsize=(10, 10))
    sns.barplot(x=feature_importances.importance, y=feature_importances.index, hue=feature_importances.importance,  dodge=False, palette='Blues_d')
    plt.gca().get_legend().set_title('Ranking')
    # invert order of the elements in the legend
    handles, labels = plt.gca().get_legend_handles_labels()
    # transform labels into %1.3f
    labels = [float(l) for l in labels]
    labels = ['%1.3f' % l for l in labels]
    plt.gca().legend(reversed(handles), reversed(labels), loc='lower right')
    plt.tight_layout(pad=1.5)
    plt.xlabel('Importance')
    plt.ylabel('Features')
    plt.title('Ranking of features')
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
    plt.figure(figsize=(10, 9.5))
    plt.plot(fpr_train, tpr_train, label='Train (area = %0.2f)' % roc_auc_train, linewidth=3)
    if test_true is not None and test_pred is not None:
        fpr_test, tpr_test, _ = roc_curve(test_true, test_pred)
        roc_auc_test = auc(fpr_test, tpr_test)
        plt.plot(fpr_test, tpr_test, label='Test (area = %0.2f)' % roc_auc_test, linewidth=3)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=20)
    plt.ylabel('True Positive Rate', fontsize=20)
    plt.title('ROC curve', fontsize=20)
    plt.legend(loc="lower right", fontsize=20)
    plt.tight_layout(pad=1.5)
    if savefig is not None:
      if savefig.endswith('.png') or savefig.endswith('.pdf'): savefig = savefig[:-4]
      for end in ['.png', '.pdf']:
        plt.savefig(savefig+end)
      plt.clf()

def DrawSigBkgHisto(train_true, train_pred, test_true=None, test_pred=None, nbins=20, savefig=None):
    plt.figure(figsize=(10, 10))
    # Normalized histograms, filled and line

    plt.hist(train_pred[train_true==1], bins=nbins, range=(0, 1), density=True, label='Signal (train)', linewidth=3, color='r', alpha=0.3)
    plt.hist(train_pred[train_true==0], bins=nbins, range=(0, 1), density=True, label='Background (train)', linewidth=3, color='b', alpha=0.3)
    if test_true is not None and test_pred is not None:
        plt.hist(test_pred[test_true==1], bins=nbins, range=(0, 1), density=True, histtype='step', label='Signal (test)', linewidth=3, color='r', alpha=0.7)
        plt.hist(test_pred[test_true==0], bins=nbins, range=(0, 1), density=True, histtype='step', label='Background (test)', linewidth=3, color='b', alpha=0.7)
    plt.xlabel('Score', fontsize=20)
    plt.ylabel('Events', fontsize=20)
    plt.title('Score distribution', fontsize=20)
    plt.legend(loc="best",  fontsize=20)
    plt.tight_layout(pad=1.5)
    if savefig is not None:
      if savefig.endswith('.png') or savefig.endswith('.pdf'): savefig = savefig[:-4]
      for end in ['.png', '.pdf']:
        plt.savefig(savefig+end)
      plt.clf()



