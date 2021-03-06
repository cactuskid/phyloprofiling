#import keras
import profiler
import os
import sys

from validation import validation_semantic_similarity

from keras.models import Sequential
from keras.layers import Dense, Activation
import pandas as pd
from utils import config_utils , hashutils
import numpy as np
# MLP for Pima Indians Dataset Serialize to JSON and HDF5
from keras.models import Sequential
from keras.layers import Dense , Dropout
from keras.models import model_from_json
from keras import optimizers , regularizers , constraints


import argparse
import glob
import time as t
import numpy
import pickle
import random
from utils import hashutils
np.random.seed(0)
random.seed(0)
import gc
###return profiler and validation obj


"""
use this script in conjuction with a labelled protein interaction datasets.
You can optimize the profile weights at diff taxonomic levels to give a
better prediction of protein interaction with the jaccard score.
A positive and negative dataset need to be generated.

"""


def load_valobjs( db , val = False , hashes=None , forest=None ):
    print('compiling' + db)
    p = profiler.Profiler(lshforestpath = forest, hashes_h5=hashes, mat_path= None , nsamples = 256, oma = True )
    print('done')
    print('loading validation')
    if val:
        if not os.path.isfile(config.datadir + 'val.pkl'):
            folder = config_utils.datadir + 'GOData/'
            val = validation_semantic_similarity.Validation_semantic_similarity( folder + 'go-basic.obo' ,
                folder + 'goframe.pkl' , folder + 'oma-go.txt' , config_utils.omadir + 'OmaServer.h5' , folder + 'termcounts.pkl' )
            with open(config.datadir + 'val.pkl' , 'wb')as valout:
                valout.write(pickle.dumps(val))
        else:

            with open(config.datadir + 'val.pkl' , 'rb')as valout:
                val = pickle.loads(valout.read())
    print( 'done')
    print('testing db')
    return p, val


def calculate_x(row):
    mat_x1 = row.mat_x
    mat_x2 = row.mat_y

    ret1 = np.zeros(mat_x1.shape)
    ret2 = np.zeros(mat_x2.shape)
    diff = mat_x1 - mat_x2
    matsum = mat_x1 + mat_x2
    ret1[np.where(diff != 0 ) ] = 1
    ret2[np.where(matsum == 2 ) ] = 1
    return list(np.hstack([ret1,ret2]))


    #generate dataframes w sem sim and profile distances

if __name__ == '__main__':



    parser = argparse.ArgumentParser()
    parser.add_argument("-hognetcsv", help="csv for training", type =str)
    parser.add_argument("-epochs", help="number of epochs to train", type=int)
    parser.add_argument("-savedir", help="save directory for model", type=str)
    parser.add_argument("-chunksize", help="num hog pairs to analyze at once", type=str)
    parser.add_argument("-overwrite", help="overwrite model", type=str)
    parser.add_argument("-name", help="name of database", type=str)
    parser.add_argument("-ntrain", help="n profiles to train on ", type=str)
    parser.add_argument("-ntest", help="n profiles to test on ", type=str)

    print(sys.argv)
    args = vars(parser.parse_args(sys.argv[1:]))

    try:
        saving_path  = self.config_utils + args['name']
    except:
        raise Exception('please specidfy db name')

    hashes_path = saving_path + 'hashes.h5'
    lshpath = saving_path + 'newlsh.pkl'
    lshforestpath = saving_path + 'newlshforest.pkl'
    mat_path = saving_path+ 'hogmat.h5'

    #load hogs dataset from paper
    #hog names need to be mapped to interactors
    if args['hognetcsv']:
        csvs = glob.glob(args['hognetcsv'] )
        print(csvs)
        df = pd.concat ( [ pd.read_csv(csv) for csv in csvs] )

    if args['epochs']:
        eps =  args['epochs']
    else:
        eps = 10
    if args['overwrite'] == 'True':
        overwrite = True
    else:
        overwrite = False


    if args['savedir']:
        savedir = args['savedir']
        if not os.path.exists(savedir):
            os.makedirs(savedir)
    else:
        savedir = './'
        # load json and create model
    if args['chunksize']:
        chunksize = args['chunksize']
    else:
        chunksize = 25

    if args['forest'] and args['hashes']:
        p = profiler.Profiler( hashes = args['hashes'], forest = args['forest'] , oma = True)
    else:
        p = profiler.Profiler( lshforestpath, hashes_path ,  oma = True)

    #shuffle
    df['HogFamA'] = df.HogA.map(hashutils.hogid2fam)
    df['HogFamB'] = df.HogB.map(hashutils.hogid2fam)

    df = df.sample(frac =1)
    msk = np.random.rand(len(df)) < 0.90

    #split
    traindf = df[msk]
    testdf = df[~msk]
    print( traindf)
    print( testdf)

    with open( config_utils.datadir + 'taxaIndex.pkl', 'rb')as taxain:
        taxaIndex = pickle.loads( taxain.read() )
    # 3 events, diff and union of matrows
    hogmat_size = 3 *  len(taxaIndex)
    if overwrite ==False & os.path.isfile(savedir + 'model.json') and  os.path.isfile(savedir + 'model.h5'):
        json_file = open(savedir + 'model.json', 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        model = model_from_json(loaded_model_json)
        # load weights into new model
        model.load_weights(savedir + "model.h5")
        print("Loaded model from disk")
    else:
        print('new model')
        layers = []
        #layers.append(Dense( 10,  , activation='relu' , use_bias=True))
        layers.append(Dense( 1, input_dim= hogmat_size , activation='elu', use_bias=False, kernel_initializer='random_uniform', kernel_constraint= constraints.NonNeg() )     )
        #layers.append( Dropout(.5 , noise_shape=None, seed=None))
        model = Sequential(layers)

    #sgd = optimizers.SGD(lr= .1, momentum=0.1, decay=0.01, nesterov=True)

    sgd = optimizers.SGD(lr= .01, momentum=0.01, decay=0.01, nesterov=True)
    #rms = optimizers.RMSprop(lr=1, rho=0.9, epsilon=None, decay=0.0)
    model.compile(loss='binary_crossentropy', optimizer=sgd, metrics=['accuracy'])

    tstart = t.time()

    X_test = np.vstack(testdf.xtrain)
    y_test = testdf.truth
    epochs = 20
    tstart= t.time()

    gendata= p.retmat_mp(traindf, nworkers = 25, chunksize=50)

    xtotal = []
    ytotal = []


    print('generate data for training')
    for i in range(int(args['ntrain'] / chunksize ) ):
        X,y next(gendata)
        xtotal.append(X)
        ytotal.append(y)
        xtotalmat = np.vstack(xtotal)
        ytotalmat = np.hstack(ytotal)
        print( xtotalmat.shape)
        print(ytotalmat.shape)
    print('training')
    metrics = model.train(x=xtotalmat  , y=ytotalmat, batch_size= 32 , epochs=1000, verbose=1 )
    print('done')

    print('generate data for testing')
    xtotal = []
    ytotal = []
    gendata= p.retmat_mp(testdf, nworkers = 25, chunksize=50)
    for i in range(int(args['ntest'] / chunksize ) ):
        X,y next(gendata)
        xtotal.append(X)
        ytotal.append(y)
        xtotalmat = np.vstack(xtotal)
        ytotalmat = np.hstack(ytotal)
        print(xtotalmat.shape)
        print(ytotalmat.shape)
    model.evaluate(x = xtotalmat , y = ytotalmat)
    print('done')
    
    print('saving')
    model_json = model.to_json()
    with open(savedir + "model_nobias.json", "w") as json_file:
        json_file.write(model_json)
    # serialize weights to HDF5
    model.save_weights(savedir + "model_nobias.h5")
    print("Saved model to disk")
    tstart = t.time()
    print('DONE')
