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
from keras.layers import Dense
from keras.models import model_from_json
import argparse
import glob
import time as t
import numpy
import pickle
import gc

###return profiler and validation obj

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

#generate dataframes w sem sim and profile distances

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-hognetcsv", help="csv for training", type =str)
    parser.add_argument("-epochs", help="number of epochs to train", type=int)
    parser.add_argument("-savedir", help="save directory for model", type=str)
    parser.add_argument("-chunksize", help="num hog pairs to analyze at once", type=str)

    parser.add_argument("-forest", help="lsh forest", type=str)
    parser.add_argument("-hashes", help="hashvals for forest", type=str)
    print(sys.argv)
    args = vars(parser.parse_args(sys.argv[1:]))
    #load hogs dataset from paper
    if args['hognetcsv']:
        csvs = glob.glob(args['hognetcsv'] )
        print(csvs)
        df = pd.concat ( [ pd.read_csv(csv) for csv in csvs] )

    if args['epochs']:
        eps =  args['epochs']
    else:
        eps = 10

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
        chunksize = 50

    if args['forest'] and args['hashes']:
        p = profiler.Profiler( hashes = args['hashes'], forest = args['forest'] , oma = True)
    else:
         p = profiler.Profiler( config_utils.datadir +'allnewlshforest.pkl', config_utils.datadir +'allhashes.h5',  oma = True)
    #shuffle
    df = df.sample(frac =1)
    msk = np.random.rand(len(df)) < 0.8

    #split
    traindf = df[msk]
    testdf = df[~msk]

    with open( config_utils.datadir + 'taxaIndex.pkl', 'rb')as taxain:
        taxaIndex = pickle.loads( taxain.read() )
    # 3 events, diff and union of matrows

    hogmat_size = 3 * 2 * len(taxaIndex)

    if os.path.isfile(savedir + 'model.json') and  os.path.isfile(savedir + 'model.h5'):
        json_file = open(savedir + 'model.json', 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(loaded_model_json)
        # load weights into new model
        loaded_model.load_weights(savedir + "model.h5")
        print("Loaded model from disk")
    else:
        model = Sequential([ Dense( 1, input_dim=2*hogmat_size , activation='linear' , use_bias=False )   ])
        model.compile(loss='mean_squared_error', optimizer='sgd', metrics=['accuracy'])

    tstart = t.time()

    testtrue = testdf.truth

    
    for batch in range(0, len(df) , chunksize ):

        slice = df.iloc[batch:batch+chunksize, :]
        interactors = list(slice.HogA) + list(slice.HogB)
        interactors = set([ hashutils.hogid2fam(h) for h in interactors])

        matrows = p.retmat_mp( interactors )
        print(matrows)
        import pdb; pdb.set_trace()

        ytrain = slice.truth
        model.train_on_batch( X_train, y_train )
        print(model.evaluate(x=X_test, y=y_test, batch_size=None, verbose=1, sample_weight=None, steps=None))
        if t.time() - tstart > 200:
            # serialize model to
            model_json = model.to_json()
            with open(savedir + "model.json", "w") as json_file:
                json_file.write(model_json)
            # serialize weights to HDF5
            model.save_weights(savedir + "model.h5")
            print("Saved model to disk")

            tstart = t.time()
