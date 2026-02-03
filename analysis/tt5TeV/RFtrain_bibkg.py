from PrepareDatasets_bibkg import *
#from NN import *
import datetime
import pandas as pd
import numpy as np
import random
import matplotlib
matplotlib.use('template')
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.decomposition import PCA
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score 
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
from pylab import rcParams
from sklearn.model_selection import train_test_split
from matplotlib.ticker import FormatStrFormatter
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import GradientBoostingClassifier
import seaborn as sns
import xgboost as xgb

from sklearn.preprocessing import label_binarize
from scipy.special import softmax

from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import HalvingRandomSearchCV
### Define path, signa, bkg, variables...

path='mva/lesscuts/'
signal = "tchannel_mva,tbarchannel_mva"
bkg = "TT_mva"#,Wjets_mva"
bkg1="Wjets_mva"

vars_train = ['2j1b_ht','2j1b_absu0eta','2j1b_sumAllPt','2j1b_leta','2j1b_j0pt','2j1b_j0eta','2j1b_u0pt','2j1b_ptjj','2j1b_mjj','2j1b_medianDRjj','2j1b_mlb', '2j1b_ptsumveclb',
	'2j1b_drlb','2j1b_ht_atlas','2j1b_mtw','2j1b_beta','2j1b_deltaeta','2j1b_topmass','2j1b_u0mass','2j1b_absleta','2j1b_topeta','2j1b_absweta','2j1b_drub','2j1b_cost2']

#['2j1b_ht','2j1b_absu0eta','2j1b_st','2j1b_sumAllPt','2j1b_leta','2j1b_j0pt','2j1b_j0eta','2j1b_u0pt','2j1b_u0eta','2j1b_ptjj','2j1b_mjj','2j1b_medianDRjj','2j1b_minDRjj','2j1b_mlb','2j1b_mt', '2j1b_ptsumveclb',
	#'2j1b_drlb','2j1b_ht_atlas','2j1b_mtw','2j1b_beta','2j1b_deltaeta','2j1b_topmass','2j1b_u0mass','2j1b_absleta','2j1b_topeta','2j1b_absweta','2j1b_drub','2j1b_cost2']  FULL LIST



### Training parameters
trainFrac = 0.85
nest=200
depth=8
split=3
leaf=3  
### Get the data
'''
signal = toList(signal); bkg = toList(bkg); vars_train = toList(vars_train)
datasetsSignal = pandas.concat([LoadColumns(path, s, vars_train, isSignal=True) for s in signal])
datasetsBkg    = pandas.concat([LoadColumns(path, b, vars_train, isSignal=False) for b in bkg], ignore_index=True)
X_train, X_test, y_train, y_test = train_test_split(datasetsBkg[vars_train], datasetsBkg['label'], test_size=0.3, random_state=4)
df_train_b=pd.concat([X_train,y_train], axis=1)
df_test_b=pd.concat([X_test,y_test], axis=1)
X_train, X_test, y_train, y_test = train_test_split(datasetsSignal[vars_train], datasetsSignal['label'], test_size=0.3, random_state=4)
df_test_g=pd.concat([X_test,y_test], axis=1)
df_train_g=pd.concat([X_train,y_train], axis=1)
#df_train_g=df_train_g.iloc[:int(df_train_g.shape[0]/2),:] 
df_train=pd.concat([df_train_b,df_train_g],axis=0)
df_test=pd.concat([df_test_b,df_test_g],axis=0)
'''
df_train, df_test = BuildDataset(path, signal, bkg,bkg1, vars_train, trainFrac, 21,nData=None)#, nData=614) this is the number of Wjets
print(df_train.sample(n=5))
### Create the model
#name='3j1b_%s_%s_minusvariablesNewMlb_train'%(nest,depth)
name='dummy_miguel'

#model=DecisionTreeClassifier(max_depth=5, min_samples_split=2, min_samples_leaf=1,class_weight='balanced')
model = RandomForestClassifier(n_estimators=nest, max_depth=depth,min_samples_split=split,min_samples_leaf=leaf,class_weight='balanced')
#model=xgb.XGBClassifier(objective='multi:softmax', random_state=21,learning_rate=0.01,n_estimators=nest,max_depth=depth )#,scale_pos_weight=2.59)
#model = GradientBoostingClassifier(n_estimators=nest, max_depth=depth)
model.n_jobs=8

#Exploring the best combination of hyperparameters:
#param_distributions= {"n_estimators":np.linspace(100,500,5).astype(int),
#                "max_depth":np.linspace(2,20,11).astype(int),
#                "min_samples_split":np.linspace(2,20,11).astype(int),
#                "min_samples_leaf":np.linspace(1,10,11).astype(int) }

#model=HalvingRandomSearchCV(model,param_distributions,verbose=1,random_state=21)

#Defining weights for xgb:
class_weights={0:1.1959,1:0.4723,2:21.4297}  #According to the formula: totalnumber/(Nclasses*Nevents in that class)
sample_weights=[class_weights[label] for label in df_train['label']] #Apply individually the weights for each event


### train!
print('training...')
model.fit(df_train[vars_train], df_train['label'])#,sample_weight=sample_weights)
print("Running on test sample. This may take a moment.")

#print('best_params',model.best_params_)
#print('best score',model.best_score_)
#print('Results',model.cv_results_)

probs = model.predict_proba(df_test[vars_train])#predict probability over test sample
probs_train = model.predict_proba(df_train[vars_train])#predict probability over test sample
#pred_y_train= model.predict(df_train[vars_train])
#pred_y= model.predict(df_test[vars_train])




print('dt test label',df_test['label'])
y_test_bin = label_binarize(df_test['label'], classes=[0, 1, 2]) #This transforms label 0,1 or 2 into [1,0,0],[0,1,0], [0,0,1] respectively
y_train_bin = label_binarize(df_train['label'], classes=[0, 1, 2]) 

print('y test bin',y_test_bin)

roc_auc=roc_auc_score(y_test_bin,probs,average='weighted')  #Esto debe de ser un tipo de media de las AUCs
print('sort of mean AUC',roc_auc)

fpr_test = dict()
tpr_test = dict()
roc_auc_test = dict()
for i in range(3):
    fpr_test[i], tpr_test[i], _ = roc_curve(y_test_bin[:, i], probs[:, i])
    roc_auc_test[i] = auc(fpr_test[i], tpr_test[i])

fpr_train = dict()
tpr_train = dict()
roc_auc_train = dict()
for i in range(3):
    fpr_train[i], tpr_train[i], _ = roc_curve(y_train_bin[:, i], probs_train[:, i])
    roc_auc_train[i] = auc(fpr_train[i], tpr_train[i])

print('roc_aucs_train',roc_auc_train)
print('roc_aucs_test',roc_auc_test)

processes=['tchannel','tt','Wjets']
# ROC

for i in range(3):
    plt.figure(figsize=(7,7))
    plt.rcParams.update({'font.size': 15}) #Larger font size
    plt.plot(fpr_test[i], tpr_test[i], color='crimson', lw=2, label='ROC curve test (area = {0:.4f})'.format(roc_auc_test[i]))
    plt.plot(fpr_train[i], tpr_train[i], color='blue', lw=1, label='ROC curve train (area = {0:.4f})'.format(roc_auc_train[i]))
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.05])
    plt.ylim([0.0, 1.05])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.title(r"$\bf{CMS}$ all vs %s"%processes[i], fontsize=20, loc='left')
    plt.savefig('mva/lesscuts/training/%s_ROCcurve_CMS.pdf'%processes[i])
    plt.savefig('mva/lesscuts/training/%s_ROCcurve_CMS.png'%processes[i])

for i in range(3):
    # probabilities
    plt.figure(figsize=(10,5))
    plt.rcParams.update({'font.size': 15}) #Larger font size
    df_test["sigprob"] = probs[:,i] #save probabilities in df
    df_train["sigprob"] = probs_train[:,i] #save probabilities in df
    
    #back = np.array(df_test["sigprob"].loc[df_test["label"]!=i].values)
    tchan = np.array(df_test["sigprob"].loc[df_test["label"]==0].values)
    tchan1 = np.array(df_train["sigprob"].loc[df_train["label"]==0].values)
    #back2 = np.array(df_train["sigprob"].loc[df_train["label"]!=1].values)
    tt = np.array(df_test["sigprob"].loc[df_test["label"]==1].values)
    tt1 = np.array(df_train["sigprob"].loc[df_train["label"]==1].values)
    #back3 = np.array(df_test["sigprob"].loc[df_test["label"]!=i].values)
    wjets = np.array(df_test["sigprob"].loc[df_test["label"]==2].values)
    wjets1 = np.array(df_train["sigprob"].loc[df_train["label"]==2].values)

    plt.hist(tchan, 20, color='blue', edgecolor='blue', lw=2, label='t-chan (test)', alpha=0.3, density=True)
    plt.hist(tchan1, 20, color='darkblue',histtype='step', lw=2, label='t-chan (train)', density=True)

    plt.hist(tt, 20, color='red', edgecolor='red', lw=2, label=r"t$\bar{t}$ (test)", alpha=0.3, density=True)
    plt.hist(tt1, 20, color='brown',histtype='step', lw=2, label=r"t$\bar{t}$ (train)",  density=True)

    plt.hist(wjets, 20, color='lightgreen', edgecolor='lightgreen', lw=2, label='Wjets (test)', alpha=0.3, density=True)
    plt.hist(wjets1, 20, color='green',histtype='step', lw=2, label='Wjets (train)',  density=True)

    plt.title(r"$\bf{CMS}$", fontsize=20, loc='left')
    plt.xlim([0, 1])
    plt.xlabel('Event probability of being classified as %s'%processes[i])
    plt.legend(loc="upper right")
    plt.grid(True)
    plt.savefig('mva/lesscuts/training/%s_probs.pdf'%processes[i])
    plt.savefig('mva/lesscuts/training/%s_probs.png'%processes[i])
    plt.show()


import pickle
pickle.dump(model, open('mva/lesscuts/training/%s_p2v2.pkl'%name, 'wb'),protocol=2) 
pickle.dump(model, open('mva/lesscuts/training/%s.pkl'%name, 'wb')) 
print('models/may/%s_p2v2.pkl'%name)

### Save the model
#from sklearn.externals import joblib
#joblib.dump(model, 'models/may/%s.joblib'%name) 
#joblib.dump(model, 'models/may/%s.pkl'%name) 
import pickle
pickle.dump(model, open('mva/lesscuts/training/%s_p2v2.pkl'%name, 'wb'),protocol=2) 
pickle.dump(model, open('mva/lesscuts/training/%s.pkl'%name, 'wb')) 
print('models/may/%s_p2v2.pkl'%name)

exit()
# ranking
feature_importances = pd.DataFrame(model.feature_importances_,index = df_test[vars_train].columns,columns=['importance']).sort_values('importance',ascending=False)                                                            
print(feature_importances)

# confusion matrix
def mostrar_resultados(y_test, pred_y):
    conf_matrix = confusion_matrix(y_test, pred_y)
    plt.figure(figsize=(5, 5))
    sns.heatmap(conf_matrix, annot=True, fmt="d");
    plt.title("Confusion matrix")
    plt.ylabel('True class')
    plt.xlabel('Predicted class')
    plt.show()
    print (classification_report(y_test, pred_y))
mostrar_resultados(df_test['label'], pred_y)
plt.savefig('mva/lesscuts/training/%s_matrix.pdf'%name)
plt.savefig('mva/lesscuts/training/%s_matrix.png'%name)


### Save the model
#from sklearn.externals import joblib
#joblib.dump(model, 'models/may/%s.joblib'%name) 
#joblib.dump(model, 'models/may/%s.pkl'%name) 
import pickle
pickle.dump(model, open('mva/lesscuts/training/%s_p2v2.pkl'%name, 'wb'),protocol=2) 
pickle.dump(model, open('mva/lesscuts/training/%s.pkl'%name, 'wb')) 
print('models/may/%s_p2v2.pkl'%name)


#sig, bkg = GetSigBkgProb(model, testdl)

exit()

#I'm refering here an alterantive way of doing this problem using DMatrix as input data instead of XGBClassifier

dtrain=xgb.DMatrix(data=df_train[vars_train],label=df_train['label'],weight=sample_weights)
dtest=xgb.DMatrix(data=df_test[vars_train],label=df_test['label'])
params = {
    'objective': 'multi:softmax',
    'num_class': 3,
    'learning_rate':0.01,
    'max_depth':depth,
    'seed': 21
}
model=xgb.train(params,dtrain,num_boost_round=nest)

#Way to define probs using the DMatrix approach
raw_margin = model.predict(dtest, output_margin=True)
probs = softmax(raw_margin, axis=1)
raw_margin_train=model.predict(dtrain,output_margin=True)
probs_train=softmax(raw_margin_train,axis=1)
