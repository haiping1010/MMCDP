import keras
import keras.backend as K
import pandas as pd
from pandas import DataFrame
from keras.models import Model, load_model
from keras.layers import Input, merge, Activation, Dropout, Dense, concatenate, Concatenate, Flatten
from keras.layers.convolutional import Convolution1D
from keras.layers.pooling import AveragePooling1D, GlobalAveragePooling1D, MaxPool1D
#from keras.layers.normalization import BatchNormalization
from keras.layers import BatchNormalization
from keras.regularizers import l2
from tensorflow.keras.optimizers import RMSprop


from tensorflow.keras.optimizers import SGD
#from keras.optimizers import RMSprop


#from keras.optimizers import SGD
from keras.callbacks import ModelCheckpoint
from keras.preprocessing.image import ImageDataGenerator

#import xgboost as xgb
from sklearn import metrics

import os
import yaml
import numpy as np
import random

random.seed(0)



def randomShuffle_name(X, Y, Z):
    random.seed(0)
    idx = [t for t in range(X.shape[0])]
    random.shuffle(idx)
    X = X[idx]
    Y = Y[idx]
    #print Z
    #print idx
    Z = Z[idx]

    print()
    print('-' * 36)
    print('dimension of X after synthesis:', X.shape)
    print('dimension of Y after synthesis', Y.shape)
    print('label after shuffle:', '\n', DataFrame(Y).head())
    print('-' * 36)
    return X, Y, Z




def randomShuffle(X, Y):
    random.seed(0)
    idx = [t for t in range(X.shape[0])]
    random.shuffle(idx)
    X = X[idx]
    Y = Y[idx]
    print()
    print('-' * 36)
    print('dimension of X after synthesis:', X.shape)
    print('dimension of Y after synthesis', Y.shape)
    print('label after shuffle:', '\n', DataFrame(Y).head())
    print('-' * 36)
    return X, Y

def synData(X_0, Y_0, X_1, Y_1,N_0,N_1, time):

    X_0_syn = X_0
    Y_0_syn = Y_0
    N_0_syn = N_0
    for i in range(time - 1):
        X_0_syn = np.vstack( (X_0_syn, X_0) )
        Y_0_syn = np.hstack( (Y_0_syn, Y_0) )
        N_0_syn =np.hstack( (N_0_syn, N_0) )

    print('dimension of generation data of X', X_0_syn.shape)
    print('dimension of generation data of Y', Y_0_syn.shape)
    print('dimension of generation data of X with label of 1', X_1.shape)
    print('dimension of generation data of Y with label of 1', Y_1.shape)

    #synthesis dataset
    X_syn = np.vstack( (X_0_syn, X_1) )
    Y_syn = np.hstack( (Y_0_syn, Y_1) )
    N_syn = np.hstack( (N_0_syn, N_1) )

    print()
    print('dimension of X after combination', X_syn.shape)
    print('dimension of Y after combination', Y_syn.shape)
    print(DataFrame(Y_syn).head())

    #shuffle data
    X_syn, Y_syn, N_syn = randomShuffle_name(X_syn, Y_syn, N_syn)

    return X_syn, Y_syn, N_syn



all_label=np.load('all_label.npy', allow_pickle=True)
all_env_vec=np.load('all_env_vec.npy', allow_pickle=True)

pos_num=len(all_env_vec)
print (pos_num)

pos=all_env_vec



test_X_pos = all_env_vec
test_X_pos = test_X_pos.astype(np.float64)



def preprocess_data_train(data_set,train_set):
    mean = np.mean(train_set)
    std = np.std(train_set)

    t = data_set

    t -= mean
    t /= std
    return t



def preprocess_data(data_set):
    #mean = np.mean(data_set, axis=0, keepdims=True)
    #std = np.std(data_set,axis=0,keepdims=True)
    #np.save('mean_std.npy', [mean, std])
    mean, std = np.load('mean_std.npy', allow_pickle=True)

    t = data_set

    t -= mean
    t /= std
    print (mean)
    print (std)
    return t



def aucJ(true_labels, predictions):
    
    fpr, tpr, thresholds = metrics.roc_curve(true_labels, predictions, pos_label=1)
    auc = metrics.auc(fpr,tpr)

    return auc

def acc(true, pred):
    
    return np.sum(true == pred) * 1.0 / len(true)

def assess(model, X, label, thre = 0.5):
    
    threshold = thre
    
    pred = model.predict(X)
    pred = pred.flatten()
    
    pred[pred > threshold] = 1
    pred[pred <= threshold] = 0
    
    auc = aucJ(label, pred)
    accuracy = acc(label, pred)
    
    print('auc: ', auc)
    print('accuracy: ', accuracy)



#train_X_old=train_X.copy()
##test_X = preprocess_data(test_X_pos )
test_X = test_X_pos
#valid_X = preprocess_data(valid_X)



def densef(X_prev, X):
    
    return concatenate([X_prev, X], axis = -1)

def dfblock(X, layer_num, dropout = 0.1, shape = 200, l2_reg = 1e-4):
    
    X_prev = Dense(shape, kernel_initializer = 'glorot_normal', activation = 'relu', kernel_regularizer = l2(l2_reg))(X)
    X_prev = BatchNormalization(axis = 1, beta_regularizer = l2(l2_reg), gamma_regularizer = l2(l2_reg))(X_prev)
    X_prev = densef(X, X_prev)
    X_prev = Dropout(dropout)(X_prev)
    
    for i in range(layer_num - 1):
        X = Dense(shape, kernel_initializer = 'glorot_normal', activation = 'relu', kernel_regularizer = l2(l2_reg))(X_prev)
        X = BatchNormalization(axis = 1, beta_regularizer = l2(l2_reg), gamma_regularizer = l2(l2_reg))(X)
        X_prev = densef(X, X_prev)
        X_prev = Dropout(dropout)(X_prev)
    
    return X_prev

def noblock(X, dropout = 0.1, shape = 800, l2_reg = 1e-4):
    
    X = Dense(shape, kernel_initializer = 'glorot_normal', activation = 'relu', kernel_regularizer = l2(l2_reg))(X)
    X = BatchNormalization(axis = 1)(X)
    X = Dropout(dropout)(X)
    
    return X

def Dense_FCNN1(input_shape, dense_layer = 3, layer_num = 3, denshape = 200, dropout = 0.1, l2_reg = 1e-4):
    
    X_input = Input(input_shape)
    
    X = dfblock(X_input, layer_num = layer_num, dropout = dropout, shape = denshape, l2_reg = l2_reg)
    X = noblock(X, dropout, input_shape[0], l2_reg)
    
    for i in range(dense_layer - 1):
        X = dfblock(X, layer_num = layer_num, dropout = dropout, shape = denshape, l2_reg = l2_reg)
        X = noblock(X, dropout, input_shape[0], l2_reg)
    
    out = Dense(21, kernel_initializer = 'glorot_normal', activation = 'sigmoid')(X)
    
    model = Model(inputs = X_input, outputs = out, name = 'FCNN1')
    
    return model

model2 = Dense_FCNN1((300, ), 1, 16, 128, 0.15, 1e-4)




model2.compile(optimizer = 'adam' , loss = "binary_crossentropy", metrics=["accuracy"])

'''
learning_rate_reduction = ReduceLROnPlateau(monitor='val_acc',
                                            patience=3,
                                            verbose=1,
                                            factor=0.5,
                                            min_lr=0.00001)
'''






model2.summary()


#checkpointer= ModelCheckpoint('model_n_{epoch:03d}.h5', verbose=1, save_weights_only=False,period=20)

#history=model2.fit(x = train_X, y = train_Y, epochs = 1500, batch_size = 64, validation_data = (valid_X,valid_Y),callbacks=[checkpointer])


model = load_model('FCdenset.h5')

#assess(model2, train_X, train_Y)
#assess(model2, valid_X, valid_Y)
#assess(model2, test_X, test_Y)


import numpy as np

aa_codes = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']




aa3to1={
   'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
   'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
   'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
   'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G'
}

aa1to3={value:key  for key,value  in aa3to1.items()}

'''
#amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']  
def get_top_3_amino_acids(predictions, amino_acids_list):
    top_3_indices = np.argsort(predictions, axis=-1)[:,-3:]
    top_3_amino_acids = np.vectorize(lambda i: amino_acids_list[i])(top_3_indices)
    return top_3_amino_acids

predictions = model.predict(test_X)
top_3_amino_acids = get_top_3_amino_acids(predictions, aa_codes)
'''

def get_amino_acids_above_threshold(predictions, amino_acids_list, threshold):
    #print predictions
    above_threshold_indices = np.where(predictions > threshold)
    #print above_threshold_indices
    above_threshold_amino_acids = np.vectorize(lambda i: amino_acids_list[i])(above_threshold_indices).flatten()
    return above_threshold_amino_acids

predictions = model.predict(test_X)

print (predictions)
dict_final={}



threshold = 0.05

for i in range( predictions.shape[0]):

   amino_acids_above_threshold = get_amino_acids_above_threshold(predictions[i], aa_codes, threshold)

   if all_label[i][0:3] in amino_acids_above_threshold:
       amino_acids_above_threshold = np.delete(amino_acids_above_threshold, np.where(amino_acids_above_threshold == all_label[i][0:3]))
   if len(amino_acids_above_threshold) != 0:
                dict_final[all_label[i][3:]] = amino_acids_above_threshold
                print (all_label[i], amino_acids_above_threshold)


#print dict_final


import pickle

with open('possible_mutations.pkl', 'wb') as f:
    pickle.dump(dict_final, f)

with open('possible_mutations.pkl', 'rb') as f:
    mutations_inf = pickle.load(f)








