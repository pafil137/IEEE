#Testing for Multi-layer Perceptron module (https://github.com/pipidog/keras_to_tensorflow/blob/master/keras_to_tensorflow.py)
from DNN import DNN_training, DNN_test, util_DNN

import pandas as pd

from keras.models import save_model

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler

##PRECONFIG####################################
#DNN - TRAINING AND TEST
isTrain = True

##DNN################################################
#DNN - CONFIGURATION
Epochs     = 20#50,
hidden_layer = [100]#[100],

activation = 'relu'
solver     = 'Adam'
metric     = 'mean_squared_error'
testPercent = 0.2#0.2/0.214/0.2062

#DATA - CONFIGURATION
configData     = 0
configDataTest = configData #None

nameData = ["ALOHA","CSMA","ALOHA+CSMA"]
data = pd.read_csv("./DATA/{0}.csv".format(nameData[configData]))

fileTrained = "./{0}.h5".format(nameData[configData]+str(-hidden_layer[0])+"x"+str(len(hidden_layer)))

wkdir = 'C:/kerasLITE/' #'C:/SAVE/kerasLITE/'
pb_filename = '{0}.pb'.format(nameData[configData]+str(-hidden_layer[0])+"x"+str(len(hidden_layer)))

##DATA###############################################
out = 'BEnh'
ent = ["h","d","R","Nb","P0","a","Pr"]
#ent = ["h","d","R","Nb", "P0", "a", "Te","Pr"]

drop = ["BEnh","Te"]

#especificacoes csma +aloha
if(configData == 2):
    data0 = pd.read_csv("./DATA/{0}.csv".format(nameData[0]))
    data1 = pd.read_csv("./DATA/{0}.csv".format(nameData[1]))
    frames = [data0,data1]
    data = pd.concat(frames)
    ent = ["h","d","R","Nb","P0","a","Pr"]
    drop = ["BEnh","Te"]
    

#REMOVER outlier
data = data.query('BEnh > -150')

print(len(data))

ydata = data[out]
Xdata = data.drop(labels = drop, axis = 1)

X, Xt, y, yt = train_test_split(Xdata, ydata, test_size=testPercent, random_state=101)

print(len(Xt))

def normalize(X_train, X_test):
    scaler     = MinMaxScaler()
    scaler.fit(X_train)
    train_keys = X_train.keys()
    test_keys  = X_test.keys()
    train_index = X_train.index
    test_index  = X_test.index
    X_train = pd.DataFrame(scaler.transform(X_train), columns=train_keys, index = train_index)
    X_test  = pd.DataFrame(scaler.transform(X_test),  columns=test_keys,  index = test_index)
    
    return X_train, X_test

if(configDataTest == configData):
    X, Xt = normalize(X, Xt)
else:
    dataT = pd.read_csv("./DATA/{0}.csv".format(nameData[configDataTest])) 
    ydataTest = dataT[out]
    XdataTest = dataT.drop(labels = ["BEnh","Pr"], axis = 1)
    #XdataTest = dataT.drop(labels = drop, axis = 1)
    X2, Xt, y2, yt = train_test_split(XdataTest, ydataTest, test_size=testPercent, random_state=101)
    X2, Xt = normalize(X, Xt)
    
##DEEP NEURAL NETWORK: FUNCTIONS####################
if(isTrain): 
    DNN, history = DNN_training.KRlt_MLPR_training(X, y, len(ent), metric, solver, hidden_layer, activation, Epochs)
    MSE = history.history[metric]
    
    #######VIEW##########
    pred=DNN.predict(Xt)
    util_DNN.litedebugDNN(list(yt), pred, MSE, "kr."+nameData[configData], '{0}'.format(nameData[configData]))
    
    #######SAVE##########  
    # HDF5 format
    save_model(DNN, fileTrained)    

else:   
    DNN = DNN_test.KRlt_MLPR_test(fileTrained)
    MSE = []
    
    pred=DNN.predict(Xt)
    util_DNN.debugDNN(list(yt), pred, X, Xt, [], nameData[configData])

'''
import math

from ModeloEnergia import _ALOHA_0  as MCE
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
import numpy as np
  
########MCE
_a = [2]
_dist = [80]
_h = [1, 2, 4, 8]

_Nb = []
_Nb.extend(range(int(math.pow(10,1)),int(math.pow(10,2)),int(math.pow(10,1))))
_Nb.extend(range(int(math.pow(10,2)),int(math.pow(10,3)),int(math.pow(10,2))))
_Nb.extend(range(int(math.pow(10,3)),int(math.pow(10,4)),int(math.pow(10,3))))
_Nb.extend(range(int(math.pow(10,4)),int(math.pow(10,5)),int(math.pow(10,4))))

_R = []
_R.extend(range(int(math.pow(10,3)),int(math.pow(10,4)),int(math.pow(10,3))))
_R.extend(range(int(math.pow(10,4)),int(math.pow(10,5)),int(math.pow(10,4))))
_R.extend(range(int(math.pow(10,5)),int(math.pow(10,6)),int(math.pow(10,5))))
_R.extend(range(int(math.pow(10,6)),int(math.pow(10,7)),int(math.pow(10,6))))

_NbS = [1000]
_RS = [1000]

BConfigMatrix, Preal = MCE.functionPowerR(_h, _dist, _R, _NbS, _a)
#BConfigMatrix, Preal = MCE.functionPowerNb(_h, _dist, _RS, _Nb, _a)

########DNN
Xmce, Ymce = util_DNN.DNNdataGeneration(BConfigMatrix)
     
Xs = StandardScaler().fit_transform(Xmce)
    
ys = np.array(Ymce)
ys = ys.flatten()                                

DNN, history = DNN_training.KRlt_MLPR_training(Xmce, ys, len(ent), metric, solver, hidden_layer, activation, Epochs) 
pred=DNN.predict(Xmce)

########COMPARATION
util_DNN.montarGraphDNNvs(BConfigMatrix, preDNN, _R, 'R(bps)')
#util_DNN.montarGraphDNNvs(BConfigMatrix, preDNN, _Nb,'Nb(bits)')
'''
