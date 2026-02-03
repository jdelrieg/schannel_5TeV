from cafea.plotter.plotter import *
import os, pandas
import numpy as np
branches = ['A_njets', 'A_nbtags', 'A_ht', 'A_st', 'A_sumAllPt', 'A_leta', 'A_j0pt', 'A_j0eta', 'A_u0pt', 'A_u0eta', 'A_ptjj', 'A_mjj', 'A_medianDRjj', 'A_minDRjj', 'A_mlb', 'A_mt', 'A_ptsumveclb', 'A_drlb']
nworkers = 64

def LoadColumns(path, sampleName, branches, isSignal=True,isWjets=False,isTT=False):
  ''' open a pkl.gz sample and retrieve columns with names "branches", that must have the same length, and returns pandas dataset '''
  path = os.path.join(path, sampleName+'.pkl.gz')
  h = loadHistos(path)
  pd = pandas.DataFrame()
  for name in branches:
    data = h[name].value
    pd.insert(0, name, data)
  if isSignal:
    pd.insert(0, 'label', np.ones(len(pd)))
  else:                                          
    pd.insert(0, 'label', np.zeros(len(pd)))
  

  return pd

def toList(t):
  ''' Transforms whatever string into a list '''
  if isinstance(t, str) and ',' in t:
    t = t.replace(' ', '').split(',')
  elif isinstance(t, str):
    t = [t]
  return t



def BuildDataset(path, signal, bkg, var=branches, trainFrac=0.8, random_state=42, nData=None, verbose=1):
    ''' Loads a dataset with labels for signal and bkg and returns train and test data samples '''
    signal = toList(signal); bkg = toList(bkg); var = toList(var);
    
    # Load datasets and apply nData per input file
    datasetsSignal = [LoadColumns(path, s, var, isSignal=True).sample(min(nData, len(LoadColumns(path, s, var, isSignal=True)))) if nData is not None else LoadColumns(path, s, var, isSignal=True) for s in signal]
    datasetsBkg = [LoadColumns(path, b, var, isSignal=False).sample(min(nData, len(LoadColumns(path, b, var, isSignal=False)))) if nData is not None else LoadColumns(path, b, var, isSignal=False) for b in bkg]
    
    # Concatenate all datasets
    datasetsSignal = pandas.concat(datasetsSignal, ignore_index=True)
    datasetsBkg = pandas.concat(datasetsBkg, ignore_index=True)
    df = pandas.concat([datasetsSignal, datasetsBkg], ignore_index=True)
    
    # Shuffle and split dataset
    train = df.sample(frac=trainFrac, random_state=random_state)
    test = df.drop(train.index)
    
    # Build info
    if verbose:
        nsig = sum(np.where(df['label'] == 1, 1, 0))
        nbkg = sum(np.where(df['label'] == 0, 1, 0))
        nsig_train = sum(np.where(train['label'] == 1, 1, 0))
        nbkg_train = sum(np.where(train['label'] == 0, 1, 0))
        
        n = len(df)
        ntrain = len(train)
        ntest = len(test)
        nvar = len(var)
        
        print(f" >> Dataset with {n} events and {nvar} columns")
        print(f" >> Contains {nsig} signal events and {nbkg} background events")
        print(f" >> The test sample has {ntest} events and the train sample has {ntrain} events")
        print(f" >> The train sample has {nsig_train} events of signal and {nbkg_train} of background")
        print(f" >> This means a signal percentage of {nsig_train / (nsig_train + nbkg_train) if (nsig_train + nbkg_train) > 0 else 0:.2f}")
    
    return train, test


def BuildDataset_balancing(path, signal, bkg, var=branches, trainFrac=0.8, random_state=42, nData=None, proportions=None, verbose=1):
    ''' Loads a dataset with labels for signal and bkg and returns train and test data samples '''
    signal = toList(signal)
    bkg = toList(bkg)
    var = toList(var)
    
    # Load data from multiple signal and background files, adding source column
    datasetsSignal = [LoadColumns(path, s, var, isSignal=True).assign(source=s) for s in signal]
    datasetsBkg = [LoadColumns(path, b, var, isSignal=False).assign(source=b) for b in bkg]
    
    df_signal = pandas.concat(datasetsSignal, ignore_index=True)
    df_bkg = pandas.concat(datasetsBkg, ignore_index=True)
    
    # Apply proportions to the signal and background datasets
    if proportions is None:
        proportions = [1/len(signal)] * len(signal) + [1/len(bkg)] * len(bkg)  # Default: equal distribution
    
    # Ensure correct length of proportions
    if len(proportions) != len(signal) + len(bkg):
        raise ValueError("Length of proportions must match the number of signal and background files.")
    
    total_events = 49066
    
    signal_samples = [df_signal[df_signal['source'] == s].sample(n=min(int(proportions[i] * total_events), len(df_signal[df_signal['source'] == s])), random_state=random_state, replace=False) for i, s in enumerate(signal)]
    bkg_samples = [df_bkg[df_bkg['source'] == b].sample(n=min(int(proportions[len(signal) + i] * total_events), len(df_bkg[df_bkg['source'] == b])), random_state=random_state, replace=False) for i, b in enumerate(bkg)]
    
    df = pandas.concat(signal_samples + bkg_samples, ignore_index=True)
    
    # Shuffle the dataset
    df = df.sample(frac=1, random_state=random_state).reset_index(drop=True)
    
    # Split into training and testing sets
    train = df.sample(frac=trainFrac, random_state=random_state)
    test = df.drop(train.index)
    
    # Display dataset information
    if verbose:
        nsig = sum(df['label'] == 1)
        nbkg = sum(df['label'] == 0)
        nsig_train = sum(train['label'] == 1)
        nbkg_train = sum(train['label'] == 0)
        n = len(df)
        ntrain = len(train)
        ntest = len(test)
        nvar = len(var)
        
        print(f" >> Dataset with {n} events and {nvar} columns")
        print(f" >> Contains {nsig} signal events and {nbkg} background events")
        print(f" >> The test sample has {ntest} events and the train sample has {ntrain} events")
        print(f" >> The train sample has {nsig_train} signal events and {nbkg_train} background events")
        print(f" >> This means a signal percentage of {nsig_train / (nsig_train + nbkg_train):.2f}")
    
    return train, test


#import torch
#from torch.utils.data import Dataset
#from torch.utils.data import DataLoader
#class HepDataset(Dataset):
#  ''' Create a torch dataloader '''
#  def __init__(self, pd):
#    self.labels = pd['label'] if 'label' in pd else np.ones(len(pd))
#    self.data   = pd.drop('label', axis=1) if 'label' in pd else pd

#  def __len__(self):
#    return len(self.data)

#  def __getitem__(self, idx):
#    dat = self.data  .iloc[idx].to_numpy(dtype=float)
#    lab = self.labels.iloc[idx]
    #lab = np.transpose([lab, np.where(lab==1, 0, 1)])
#    return dat, lab

#def PrepareData(path, signal, bkg, var=branches, trainFrac=0.8, batch_size=64, nData=None, verbose=1):
#  dtrain, dtest = BuildDataset(path, signal, bkg, var, trainFrac, nData, verbose)
#  train    = HepDataset(dtrain)
#  train_dl = DataLoader(train, batch_size=batch_size, shuffle=True, num_workers=nworkers)
#  test     = HepDataset(dtest)
#  test_dl  = DataLoader(test, batch_size=batch_size, shuffle=True, num_workers=nworkers)
#  return train_dl, test_dl
