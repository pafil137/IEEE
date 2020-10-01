#TRAINING for Multi-layer Perceptron module (sklearn.neural_network)
#RADIO COGNITIVO

from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_squared_error

import tensorflow as tf

from tensorflow.python.keras.models import Sequential
from tensorflow.python.keras.layers import Dense

from keras.models import Sequential as KSequential
from keras.layers import Dense as KDense

##PROCESSAMENTO##############################
#DEEP
def SK_MLPR_training(X, y, minMSE, momentumAll, solver, hidden_layer, activation, learning_rate, qtd_batch):
    print("####CRIANDO MLPR####")
    
    for momentum in momentumAll:
        
        mlp = MLPRegressor(solver=solver, 
                           hidden_layer_sizes = hidden_layer,
                           activation=activation,
                           learning_rate_init=learning_rate, 
                           random_state=1,
                           batch_size=qtd_batch, 
                           momentum=momentum)
        
        #TREINO
        mse = []
        bmse = 99999
        i = 0
        
        #for i in range(Epochs):
        while bmse > minMSE:
            i+=1
            mlp.partial_fit(X, y)
            bmse = mean_squared_error(y, mlp.predict(X))
            mse.append(bmse)
            print(i," M.S.E:",bmse)
        
        return mlp, mse

def TF_MLPR_training(training_input_fn, linear_features, solver, hidden_layer, activation, Epochs, fileTrained):
    print("####CRIANDO MLPR####")
    sess = tf.Session()
    DNN  = tf.estimator.DNNRegressor(hidden_units = hidden_layer,
                                feature_columns   = linear_features,
                                activation_fn     = activation,
                                optimizer         = solver,
                                model_dir         = fileTrained)     
      
    DNN.train(input_fn = training_input_fn, steps = Epochs)
    
    sess.close()
    
    return DNN

def KR_MLPR_training(X, y, input_dim, metric, solver, hidden_layer, activation, Epochs, qtd_batch):
    
    DNN = Sequential()
    DNN.add(Dense(hidden_layer, input_dim=input_dim, kernel_initializer='normal', activation=activation))
    #DNN.add(Dense(hidden_layer, activation=activation))
    #DNN.add(Dense(hidden_layer, activation=activation))
    DNN.add(Dense(1, activation='linear'))
    DNN.summary()
    
    DNN.compile(loss=metric, optimizer=solver, metrics=[metric,'mae'])
    history = DNN.fit(X, y, epochs=Epochs, batch_size=qtd_batch,  verbose=1)
    
    return DNN, history
    
def KRlt_MLPR_training(X, y, input_dim, metric, solver, hidden_layer, activation, Epochs):
    
    DNN = KSequential()
    
    DNN.add(KDense(units = input_dim, input_shape=(input_dim,)))
    
    for layer in hidden_layer:
        DNN.add(KDense(units = layer, activation =activation))
    
    DNN.add(KDense(units = 1, activation ='linear'))
    
    DNN.compile(loss=metric, optimizer=solver, metrics=[metric])
    history = DNN.fit(x = X, y=y, epochs = Epochs) 
    
    return DNN, history
  
def KRlt_MLPR_trainingREV(X, y, input_dim, output_dim, metric, solver, hidden_layer, activation, Epochs):
    
    DNN = KSequential()
    DNN.add(KDense(units = hidden_layer, input_shape=(input_dim,), activation =activation))
    DNN.add(KDense(units = hidden_layer, activation =activation))
    DNN.add(KDense(units = hidden_layer, activation =activation))
    DNN.add(KDense(units = output_dim, activation ="linear"))
    
    DNN.compile(loss=metric, optimizer=solver, metrics=[metric])
    history = DNN.fit(x = X, y=y, epochs = Epochs) 
    
    return DNN, history

