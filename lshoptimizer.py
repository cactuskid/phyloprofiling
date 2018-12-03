from lshbuilder import LSHBuilder
import profiler
from validation import validation_semantic_similarity
from utils import config_utils
from functools import partial
import pickle
import functools
import itertools
from bayes_opt import BayesianOptimization
from utils import hashutils
import numpy as np
import random

#set rand seed
np.random.seed(0)
random.seed(0)
from tables import *

def profiling_error( db , taxfilter, tax_mask, lossweight , presenceweight, dupweight, loss_lambda , presence_lambda , dupl_lamba,  hoglist ):
    print('compiling' + db)
    #record param settings
    #compile lsh
    parastr = 'lw'+str(lossweight)+ 'pw'+str(presenceweight)+ 'dw'+ str(dupweight)+ 'll'+str(loss_lambda)+ 'pl'+ str(presenceweight) +'dl' + str(dupl_lamba)
    #print(parastr)
    startdict={'presence':presenceweight, 'loss':lossweight, 'dup':dupweight}
    lambdadict={'presence':presence_lambda, 'loss':loss_lambda, 'dup':dupl_lamba}

    with open_file(config_utils.omadir + 'OmaServer.h5', mode="r") as h5_oma:
        lsh_builder = LSHBuilder(h5_oma, saving_folder= config_utils.datadir , saving_name=db, numperm = 256,
        treeweights= None , taxfilter = taxfilter, taxmask= tax_mask , lambdadict= lambdadict, start= startdict)
        hashes, forest , mat = lsh_builder.run_pipeline()
        #hashes, forest, lshpath =lsh_builder.run_pipeline()
    print( 'done compiling')
    print('query DB and calculate error')
    print('load profiler')
    p = profiler.Profiler(lshforestpath = forest, hashes_h5=hashes, mat_path= None )
    print('done')
    print('loading validation')

    folder = config_utils.datadir + 'GOData/'
    val = val = validation_semantic_similarity.Validation_semantic_similarity( folder + 'go-basic.obo' ,
        folder + 'goframe.pkl' , folder + 'oma-go.txt' , config_utils.omadir + 'OmaServer.h5' , folder + 'termcounts.pkl' )

    print( 'done')
    print('testing db')

    if not hoglist:
        #sample random hogs
        hoglist = list(np.random.randint(0, high=610000, size=5000, dtype='l'))
        hoglist = [ hashutils.fam2hogid(s)  for s in hoglist]

    scores = {}
    for i,hog in enumerate(hoglist):
        print(hog)
        res = p.hog_query( hog_query= hog , k = 20)
        res = [ hashutils.fam2hogid(r) for r in res]
        scores.update( { combo: {'query_num':i, 'hog_sem_sim': val.semantic_similarity_score(combo[0], combo[1])
        , 'hog_resnik_sim' : p.hog_v_hog(combo[0], combo[1])
        }  for combo in itertools.combinations(res,2) } )
    
    resdf = pd.DataFrame.from_dict( scores, orient = 'index')
    resdf.to_csv( config_utils.datadir + 'resdf_' +  parastr + '.csv')
    #take positive information values
    errorval = resdf[resdf.hog_sem_sim >0].hog_sem_sim.mean()
    print(errorval)
    print('done')
    print(errorval)
    return errorval


class Observer_saver():
    def __init__(self, savepath):
        pass
    def update(self, event):

        if event is Events.INIT_DONE:
            print("Initialization completed")
        elif event is Events.FIT_STEP_DONE:
            print("Optimization step finished, current max: ", instance.res['max'])
        elif event is Events.FIT_DONE:
            print("Optimization finished, maximum value at: ", instance.res['max'])

if __name__ == '__main__':

    dbdict = {
    #'plants': { 'taxfilter': None , 'taxmask': 33090 }
    'all': { 'taxfilter': None , 'taxmask': None }
    #'archaea':{ 'taxfilter': None , 'taxmask': 2157 }
    #'bacteria':{ 'taxfilter': None , 'taxmask': 2 }
    #'eukarya':{ 'taxfilter': None , 'taxmask': 2759 }
    #'protists':{ 'taxfilter': [2 , 2157 , 33090 , 4751, 33208] , 'taxmask':None }
    #'fungi':{ 'taxfilter': None , 'taxmask': 4751 }
    #'metazoa':{ 'taxfilter': None , 'taxmask': 33208 }
    }

    for db in dbdict:
        print(db)
        #compile validation HOG list

        num_rounds = 3000
        random_state = 2016
        num_iter = 25
        init_points = 5

        observer = Observer_saver('./')

        #give the validation hog list and taxfilter and taxmask parameters

        error = functools.partial( profiling_error , db=db , taxfilter = dbdict[db]['taxfilter'], tax_mask = dbdict[db]['taxmask'],  hoglist =None)

        #get error for the first point with all weights at 1
        bo = BayesianOptimization( error ,  {'lossweight': (0, 1),
                                                        'presenceweight': (0, 1),
                                                        'dupweight':(0,1),
                                                        'loss_lambda':(-1,1),
                                                        'presence_lambda':(-1,1),
                                                        'dupl_lamba':(-1,1)
                                                        })
        import pdb; pdb.set_trace()
        bo.maximize(init_points=5, n_iter=15, kappa=2)

        #save the friggin result
        with open( './bayesopt.pkl', mode='wb', buffering=None) as bayesout:
            bayesout.write(pickle.dumps(bo, -1))

        print(bo.res['max'])

        # Making changes to the gaussian process can impact the algorithm
        # dramatically.
        gp_params = {'kernel': None,
                     'alpha': 1e-5}

        # Run it again with different acquisition function
        bo.maximize(n_iter=5, acq='ei', **gp_params)

        # Finally, we take a look at the final results.
        print(bo.res['max'])
        print(bo.res['all'])
