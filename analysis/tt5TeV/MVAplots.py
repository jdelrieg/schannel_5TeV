from cafea.plotter.DataGraphs import *
import pickle as pkl
from config import *



pathmodel = '/mnt_pool/c3_users/user/jriego/tchannel5TeV/splitting_tchan_tbar/mva/all_stats_fix/training_feb/b10_unbalanced.pkl' 

random_state = 21 # 42
vars2j1b = ['2j1b_absu0eta','2j1b_deltaeta','2j1b_topmass','2j1b_mt','2j1b_medianDRjj','2j1b_mlb','2j1b_mjj','2j1b_mub','2j1b_sumAllPt','2j1b_ht_atlas','2j1b_DptWub','2j1b_DphiWub','2j1b_Detalu']#,'2j1b_j0eta','2j1b_u0pt','2j1b_u0eta','2j1b_ptjj','2j1b_mjj','2j1b_medianDRjj','2j1b_minDRjj','2j1b_mlb','2j1b_mt', '2j1b_ptsumveclb']


#baseweb='/nfs/fanae/user/jriego/www/public/tt5TeV/'
outpath = 'mva/all_stats_fix/training_feb/furtherplots/'

path = '/mnt_pool/c3_users/user/jriego/tchannel5TeV/splitting_tchan_tbar/mva/all_stats_fix/'
signal = ["tchannel_mva,tbarchannel_mva"]
bkg = ["TT_mva,Wjets_mva"]#"W0JetsToLNu, W1JetsToLNu, W2JetsToLNu, W3JetsToLNu"#WJetsToLNu, 

datadict = {'signal': signal, 'bkg': bkg}

def CheckMaxMin(modelpath, df, vars):
  model = pkl.load(open(modelpath, 'rb'))
  pred = model.predict_proba(df[vars])[:,1]
  sig = pred[df['label']==1]
  bkg = pred[df['label']==0]
  print('Min and max for signal     = (%1.2f, %1.2f)'%(min(sig), max(sig)))
  print('Min and max for background = (%1.2f, %1.2f)'%(min(bkg), max(bkg)))
def DrawPlotsModel(modelpath, vars, lev):
    outpathlev = f'{outpath}{lev}/'
    print(outpathlev)
    
    if not os.path.exists(outpathlev): os.makedirs(outpathlev)
    ### Get the dataframe and draw some characteristic plots
    df = BuildPandasDF(path, datadict, vars, even=False)    #Esto iguala entre fondo y señal (si even=True)
    CheckMaxMin(modelpath, df, vars)
    DrawHistos(df, vars2j1b, savefig=f'{outpathlev}histos_{lev}.png')
    DrawBoxPlots(df, vars, savefig=f'{outpathlev}boxplots_{lev}.png')
    DrawCorrelationMatrix(df, vars, savefig=f'{outpathlev}correlation_{lev}.png')


    # Load the model
    model = pkl.load(open(modelpath, 'rb'))

    ### Draw ranking
    rank = DrawRanking(model, df, vars, savefig=f'{outpathlev}ranking_{lev}.png')
    #DrawPairPlots(df, list(rank.index)[:5], savefig=f'{outpathlev}pairplots_{lev}.png')

    # Confusion matrix
    y_true = df['label'].values
    y_pred = model.predict(df[vars])
    ConfusionMatrix(y_true, y_pred, savefig=f'{outpathlev}confusion_{lev}.png')

    # Histogram of probabilities for signal and background and ROC curve
    df_train, df_test = BuildPandasDF(path, datadict, vars, even=True, train_test_split=0.85, random_state=random_state)
    train_true = df_train['label'].values
    train_pred = model.predict_proba(df_train[vars])[:,1]
    test_true = df_test['label'].values
    test_pred = model.predict_proba(df_test[vars])[:,1]

    print(len(df_train),len(df_test))
    #print(df.head(5),'\n',model.predict_proba(df[vars])[:,1][:5])

    DrawROC(train_true, train_pred, test_true, test_pred, savefig=f'{outpathlev}roc_{lev}.png')
    DrawSigBkgHisto(train_true, train_pred, test_true, test_pred, savefig=f'{outpathlev}sigbkg_{lev}.png')
    
    prediccion=model.predict_proba(df[vars])[:,1]
    corte_barbara=(prediccion>0.3)&(prediccion<0.55)
    df_barbara=df[corte_barbara]
    df_barbara_anti=df[~corte_barbara]
    
    DrawHistos_double(df_barbara_anti,df_barbara, vars2j1b, savefig=f'{outpathlev}histos_double_{lev}.png')
    #DrawHistos(df_barbara_anti, vars2j1b, savefig=f'{outpathlev}histos_sides_{lev}.png')

if __name__ == '__main__':
  # Training plots
  DrawPlotsModel(pathmodel, vars2j1b, '2j1b')
