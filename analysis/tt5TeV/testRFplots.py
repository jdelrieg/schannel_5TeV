from cafea.plotter.DataGraphs import *
import pickle as pkl

#pathmodels = '/nfs/fanae/user/andreatf/cafea/cafea/analysis/tt5TeV/models/november/'
pathmodels = '/mnt_pool/c3_users/user/andreatf/cafea/cafea/analysis/tt5TeV/models/january/'
random_state = 50 # 42



#model3j1b = 'rf3j1b_5000_6_allvariablesNewMlb.pkl'  # 'rf3j1b_200_4_allvariablesNewMlb_p2v2.pkl' #rf3j1b_200_7_allvariables_p2v2.pkl
model3j1b = '3j1b_500_6_minusvariablesNewMlb.pkl'
#vars3j1b = ['3j1b_ht', '3j1b_st', '3j1b_sumAllPt', '3j1b_j0pt', '3j1b_u0pt', '3j1b_ptjj', '3j1b_mjj', '3j1b_medianDRjj', '3j1b_minDRjj', '3j1b_mlb', '3j1b_ptsumveclb', '3j1b_drlb', '3j1b_druumedian', '3j1b_muu', '3j1b_ptuu']
vars3j1b = ['3j1b_ht','3j1b_j0pt', '3j1b_mjj', '3j1b_medianDRjj', '3j1b_mlb', '3j1b_drlb', '3j1b_druumedian', '3j1b_muu']

#model3j2b = 'rf3j2b_250_4_allvariablesNewMlb_p2v2.pkl'
#vars3j2b = ['3j2b_ht', '3j2b_st', '3j2b_sumAllPt', '3j2b_j0pt', '3j2b_j0eta', '3j2b_ptjj', '3j2b_mjj', '3j2b_medianDRjj', '3j2b_minDRjj', '3j2b_mlb', '3j2b_ptsumveclb', '3j2b_drlb']


outpath = '/nfs/fanae/user/juanr/www/public/tt5TeV/ljets/11jan2023/RF/'

path = '/mnt_pool/c3_users/user/juanr/cafea/histos5TeV/28nov2022/forTraining/'
signal = ["TTPS"]
bkg = ["WJetsToLNu", "W0JetsToLNu", "W1JetsToLNu", "W2JetsToLNu", "W3JetsToLNu"]

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
    if not os.path.exists(outpathlev): os.mkdir(outpathlev)
    ### Get the dataframe and draw some characteristic plots
    df = BuildPandasDF(path, datadict, vars, even=True)
    CheckMaxMin(modelpath, df, vars)
    DrawHistos(df, vars, savefig=f'{outpathlev}histos_{lev}.png')
    DrawBoxPlots(df, vars, savefig=f'{outpathlev}boxplots_{lev}.png')
    DrawCorrelationMatrix(df, vars, savefig=f'{outpathlev}correlation_{lev}.png')


    # Load the model
    model = pkl.load(open(modelpath, 'rb'))

    ### Draw ranking
    rank = DrawRanking(model, df, vars, savefig=f'{outpathlev}ranking_{lev}.png')
    DrawPairPlots(df, list(rank.index)[:5], savefig=f'{outpathlev}pairplots_{lev}.png')

    # Confusion matrix
    y_true = df['label'].values
    y_pred = model.predict(df[vars])
    ConfusionMatrix(y_true, y_pred, savefig=f'{outpathlev}confusion_{lev}.png')

    # Histogram of probabilities for signal and background and ROC curve
    df_train, df_test = BuildPandasDF(path, datadict, vars, even=True, train_test_split=0.7, random_state=random_state)
    train_true = df_train['label'].values
    train_pred = model.predict_proba(df_train[vars])[:,1]
    test_true = df_test['label'].values
    test_pred = model.predict_proba(df_test[vars])[:,1]

    DrawROC(train_true, train_pred, test_true, test_pred, savefig=f'{outpathlev}roc_{lev}.png')
    DrawSigBkgHisto(train_true, train_pred, test_true, test_pred, savefig=f'{outpathlev}sigbkg_{lev}.png')

DrawPlotsModel(pathmodels+model3j1b, vars3j1b, '3j1b')
#DrawPlotsModel(pathmodels+model3j2b, vars3j2b, '3j2b')
