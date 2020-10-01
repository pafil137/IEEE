#TESTING for Multi-layer Perceptron module (sklearn.neural_network)

#MLP
from sklearn.neural_network import MLPRegressor
import _pickle as cPickle

import tensorflow as tf
from tensorflow.python.keras.models import load_model

from keras.models import load_model as Kload

#CARREGAR MLPR SALVO EM ARQUIVO        
###SKLEARN##############################################
def SK_MLPR_test(Xt, DNN_loaded):
    return DNN_loaded.predict(Xt)


###TENSORFLOW###########################################
def TF_MLPR_test(linear_features, solver, hidden_layer, activation, Epochs, fileTrained):
    sess = tf.Session()
    DNN = tf.estimator.DNNRegressor(hidden_units    = hidden_layer,
                                feature_columns = linear_features,
                                activation_fn   = activation,
                                optimizer       = solver,
                                model_dir       = fileTrained) 
    
    saver = tf.train.import_meta_graph(fileTrained+"model.ckpt-{0}.meta".format(Epochs))
    saver.restore(sess,tf.train.latest_checkpoint(fileTrained))
    
    sess.close()
    
    return DNN

###KERAS################################################
def KR_MLPR_test(fileTrained):
    return load_model(fileTrained)

###KERAS-lite###########################################
def KRlt_MLPR_test(fileTrained):
    return Kload(fileTrained)