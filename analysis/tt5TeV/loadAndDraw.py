import pickle
from DrawModel import *
from NN import EvaluateModelForDataset

model_name = 'models/model_04Jul22_07h04m.pkl'

### Open the model
with open(model_name, 'rb') as f:
  model = pickle.load(f)
  
path = '../histos/tt5TeV/forNN/'
signal = "TT, TTPS"
bkg = "WJetsToLNu, W0JetsToLNu, W1JetsToLNu, W2JetsToLNu, W3JetsToLNu"
var = ['A_ht', 'A_sumAllPt', 'A_leta', 'A_j0pt', 'A_mjj', 'A_medianDRjj', 'A_drlb']

trainFrac = 0.8
batch_size = 2048
learning_rate = 0.001
epochs = 200

traindl, testdl = PrepareData(path, signal, bkg, var, trainFrac, batch_size, nData=-1)
print(model)
prob_test , lab_test  = EvaluateModelForDataset(model, testdl )
prob_train, lab_train = EvaluateModelForDataset(model, traindl)

PlotROC(prob_train, lab_train, prob_test, lab_test)
PlotHisto(prob_train, lab_train)
