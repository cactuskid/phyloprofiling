import keras
import profiler
import os
from validation import Validation_semantic_similarity

from keras.models import Sequential
from keras.layers import Dense, Activation

# MLP for Pima Indians Dataset Serialize to JSON and HDF5
from keras.models import Sequential
from keras.layers import Dense
from keras.models import model_from_json
import time as t
import numpy
import os

###return profiler and validation obj
def load_valobjs( db , val = False , hashes=None , forest=None ):
    print('compiling' + db)
    p = profiler.Profiler(lshforestpath = forest, hashes_h5=hashes, mat_path= None , nsamples = 256)
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
    print(sys.argv)
    args = vars(parser.parse_args(sys.argv[1:]))
    #load hogs dataset from paper
    if args['hognetcsv']:
        csv =  args['hognetcsv']
        df = pd.read_csv(csv)
    if args['epochs']:
        eps =  args['epochs']
    else:
        eps = 10


    if arsgs['savedir']:
        savedir = args['savedir']
        if not os.path.exists(savedir):
            os.makedirs(savedir)
    else:
        savedir = './'
        # load json and create model

    if os.path.isfile(savedir + 'model.json') and  os.path.isfile(savedir + 'model.h5')
        json_file = open(savedir + 'model.json', 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(loaded_model_json)
        # load weights into new model
        loaded_model.load_weights(savedir + "model.h5")
        print("Loaded model from disk")
    else:
        model = Sequential([ Dense( 1, input_dim=2*hogmat_size ) , Activation('relu')  use_bias=True  ])
        model.compile(loss='mean_squared_error', optimizer='sgd', metrics=['accuracy'])
    tstart = t.time()

    interactors =  list(df.HogA) + list(df.HogB)
    iteractors = set([ hashutils.hog2famid(h) for h in interactors])
    matrows = p.retmat_mp( interactors )

    df['HogA_row'] = df.HogA.map()

    for batch in generator:

        x = np.vstack(df.rows)
        y = np.asarray(df.truth)

        X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.33, random_state=0)
        model.train_on_batch( X_train, y_train )
        print(model.evaluate(x=X_test, y=y_test, batch_size=None, verbose=1, sample_weight=None, steps=None))
        if t.time() - ttart > 200:
            # serialize model to JSON
            model_json = model.to_json()
            with open(savedir + "model.json", "w") as json_file:
                json_file.write(model_json)
            # serialize weights to HDF5
            model.save_weights(savedir + "model.h5")
            print("Saved model to disk")
