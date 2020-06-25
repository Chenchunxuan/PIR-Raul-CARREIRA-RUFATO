# -*- coding: utf-8 -*-
"""
Created on Fri May 29 09:25:30 2020

@author: RaulCarreira

    This routine creates and trains a neural network using an input data 
    (external txt file).
    
    1) To change the input file path, change modify line 73
    2) To see the Kfold validation score, uncomment from line 134 to line 137
    3) To change the model, change function build_Model()
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import tensorflow.keras
from sklearn.model_selection import train_test_split
import pandas as pd
from keras import layers
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras
import keras.layers.core
from sklearn.preprocessing import StandardScaler
# import numpy as np
from sklearn.model_selection import KFold,cross_val_score
from keras.wrappers.scikit_learn import KerasRegressor
from sklearn import preprocessing
from keras.layers import LeakyReLU
# import pandas as pd
from keras.optimizers import adam
from keras.layers import Layer
from keras import backend as K
import GenerateAirfoil as pl
from sklearn import linear_model
from tensorflow.keras.callbacks import EarlyStopping
from scipy.stats import pearsonr
import numpy as np
from sklearn.linear_model import LinearRegression

class RBFLayer(Layer):
    # This class creates a layer activation function to use: Radial basis function
    def __init__(self, units, gamma, **kwargs):
        super(RBFLayer, self).__init__(**kwargs)
        self.units = units
        self.gamma = K.cast_to_floatx(gamma)

    def build(self, input_shape):
        self.mu = self.add_weight(name='mu',
                                  shape=(int(input_shape[1]), self.units),
                                  initializer='uniform',
                                  trainable=True)
        super(RBFLayer, self).build(input_shape)

    def call(self, inputs):
        diff = K.expand_dims(inputs) - self.mu
        l2 = K.sum(K.pow(diff, 2), axis=1)
        res = K.exp(-1 * self.gamma * l2)
        return res

    def compute_output_shape(self, input_shape):
        return (input_shape[0], self.units)

print(tf.__version__)


def remove_outlier(df_in, col_name):
    # df_in: .csv file desired to change
    # col_name:  name of the column to change
    
    """ This function allows to limit the extreme values ​​in a Gaussian distribution
        
    Exemple): If you want to cut the gaussian distribution in a superior limit
    of 80% and an inferior limit of 20%:
        
        q1 = df_in[col_name].quantile(0.2)
        q3 = df_in[col_name].quantile(0.8)
    
    """
    q1 = df_in[col_name].quantile(0.01)
    q3 = df_in[col_name].quantile(0.99)
    iqr = q3-q1 #Interquartile range
    fence_low  = q1-1.5*iqr
    fence_high = q3+1.5*iqr
    df_out = df_in.loc[(df_in[col_name] > fence_low) & (df_in[col_name] < fence_high)]
    return df_out


archive2 = open('Training.txt', 'r')            # load dataset

read_file = pd.read_csv (r'Training.txt')       # transform into a .csv file
read_file.to_csv (r'Training.csv', index=None)
white = pd.read_csv('Training.csv',sep=';')     # cuting the .csv file in ";"

white.isnull().any()                            # excluding bad data
white = white.dropna()                          # excluding bad data

for (columnName) in white.iteritems():
    # Performing the cut on the Gaussian curve
    df_out = remove_outlier(white, columnName[0])

white = df_out



min_max_scaler = preprocessing.MinMaxScaler()
new = min_max_scaler.fit_transform(white)

# Specify the inputs
X = white.iloc[:,0:31]

#Specify the outputs 
y = white.iloc[:,31:46]


# Generating train data and test data randomly (test_size in %)
X_train, X_test, y_train_first, y_test_first = train_test_split(X, y, test_size=0.1, random_state=42)


normed_train_data = X_train.values
normed_test_data  = X_test.values
y_train  =  y_train_first.values
y_test = y_test_first.values



# min_max_scaler = preprocessing.MinMaxScaler()
# normed_train_data = min_max_scaler.fit_transform(X_train)

# min_max_scaler2 = preprocessing.MinMaxScaler()
# y_train = min_max_scaler2.fit_transform(y_train_first)

# min_max_scaler = preprocessing.MinMaxScaler()
# normed_test_data = min_max_scaler.fit_transform(X_test)

# min_max_scaler2 = preprocessing.MinMaxScaler()
# y_test = min_max_scaler2.fit_transform(y_test_first)
        
# sc_X = StandardScaler()
# normed_train_data = sc_X.fit_transform(X_train)

# sc_X2 = StandardScaler()
# y_train = sc_X2.fit_transform(y_train_first)

#print(np.any(np.isnan(normed_test_data)))
#print(np.any(np.isnan(y_test)))
# LeakyReLU(alpha=0.01)

def build_model():
    "Function that builds the model"
    model = keras.Sequential([
        layers.Dense(500, activation= LeakyReLU(alpha=0.01), input_shape=[len(X_train.keys())]),
        layers.Dropout(0.5),
        # RBFLayer(100,0.9),
        layers.Dense(100, activation = LeakyReLU(alpha=0.01)), 
        layers.Dropout(0.5),
        layers.Dense(50, activation  = LeakyReLU(alpha=0.01)),
        layers.Dropout(0.5),
        layers.Dense(15)
    ])

  # optimizer = Adam(learning_rate=0.001,beta_1=0.999, beta_2=0.999,epsilon=1e-07,amsgrad=False,name="Adam")

    optimizer = adam(learning_rate=0.001, beta_1=0.99, beta_2=0.999,amsgrad=False)#,clipnorm=0.001
    model.compile(loss='mse', optimizer=optimizer, metrics=['mae', 'mse'])
    return model


# "To do a K-fold validation uncomment from line 153-156"
# estimator = KerasRegressor(build_fn=build_model, epochs=100, batch_size=100, verbose=1)
# kfold = KFold(n_splits=10)
# results = cross_val_score(estimator, normed_train_data, y_train, cv=kfold)
# print("Baseline: %.8f (%.8f) MSE" % (results.mean(), results.std()))


model = build_model()
model.summary()
# monitor = EarlyStopping(monitor='val_loss', mode='min')
# monitor = EarlyStopping(monitor='val_loss', min_delta=0.0001, 
#                         patience=300, verbose=1, mode='min', 
#                         restore_best_weights=True)

history = model.fit(normed_train_data, y_train,epochs=  300, batch_size=50,
                    validation_data = (normed_test_data,y_test), verbose=2) #,callbacks=[monitor]

hist = pd.DataFrame(history.history)
hist['epoch'] = history.epoch
hist.tail()

# model.save("model.h5")
# print("Saved model to disk")


def plot_history(history):
    "Funtion that plots the errors history"
    hist = pd.DataFrame(history.history)
    hist['epoch'] = history.epoch
      
    plt.figure()
    plt.xlabel('Epoch')
    plt.ylabel('Mean Abs Error [MPG]')
    plt.yscale('log')
    plt.plot(hist['epoch'], hist['mae'],label='Train Error')
    plt.plot(hist['epoch'], hist['val_mae'],label = 'Test Error')
    # plt.ylim([0,5])
    plt.legend()
    plt.grid()
      
    plt.figure()
    plt.xlabel('Epoch')
    plt.ylabel('Mean Square Error [$MPG^2$]')
    plt.yscale('log')
    plt.plot(hist['epoch'], hist['mse'],label='Train Error')
    plt.plot(hist['epoch'], hist['val_mse'],label = 'Test Error')
    # plt.ylim([0,20])
    plt.legend()
    plt.grid()
    plt.show()

def chart_regression(pred, y, sort = True):
    "Doing a chart_regressiont to see the outliers "
    t = pd.DataFrame({'pred': pred, 'y': y.flatten()})
    if sort:
        t.sort_values(by=['y'], inplace=True)
    plt.plot(t['y'].tolist(), label='expected')
    plt.plot(t['pred'].tolist(), label='prediction')
    plt.ylabel('output')
    plt.legend()
    plt.show()
    

plot_history(history)

"Plotting predictions and fitted train data"
plt.figure()
test_predictions = model.predict(normed_train_data[170:250]).flatten()
plt.scatter(y_train[170:250], test_predictions, s=50, facecolors='none', edgecolors='r')
plt.xlabel('True Values [MPG]')
plt.ylabel('Predictions [MPG]')
plt.axis('equal')
plt.axis('square')
plt.xlim([-0.5,1.5])
plt.grid()
plt.ylim([-0.5,1.5])
_ = plt.plot([-100, 100], [-100, 100],'b')
plt.show()

plt.figure()
test_predictions = model.predict(normed_test_data).flatten()
y_test_vector = np.zeros(len(test_predictions))
c=0
for i in range(len(y_test[:,0])):
    for j in range(len(y_test[0,:])):
        y_test_vector[j+i*len(y_test[0,:])] = y_test[i,j]
corr, _ = pearsonr(y_test_vector,test_predictions)

lm = linear_model.LinearRegression() # Linear regression to take Clalpha
model2 = lm.fit(y_test_vector.reshape(-1, 1),test_predictions)
alpha = model2.coef_
alpha0 = model2.intercept_
er = np.linspace(-10,10,200)
curv = alpha0+alpha*er


plt.scatter(y_test, test_predictions, s=30, facecolors='none', edgecolors='r')
plt.plot(er, curv,'y')
plt.xlabel('True Values [MPG]')
plt.ylabel('Predictions [MPG]')
plt.axis('equal')
plt.axis('square')
plt.xlim([-0.5,1.5])
plt.grid()
plt.ylim([-0.5,1.5])
_ = plt.plot([-100, 100], [-100, 100],'b')

# # Plot the chart
# plt.figure()
# chart_regression(test_predictions,y_test)


        




def Find_Airfoil(desired,n):
    # Desired vector(must have the same size as trainning data)
    "After train your N.N., this functions allows to use it"
    Vector = model.predict(desired).flatten()
    print(Vector)
    X,Y = pl.BuildAirfoil(Vector,n)
    return X,Y

