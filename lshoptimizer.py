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
import multiprocessing as mp
import time

#set rand seed
np.random.seed(0)
random.seed(0)
from tables import *

def profiling_error( db , taxfilter, tax_mask, lossweight , presenceweight, dupweight, loss_lambda , presence_lambda , dupl_lamba,  hoglist , val = None, compile = True):
    print('compiling' + db)
    #record param settings
    #compile lsh
    parastr = 'lw'+str(lossweight)+ 'pw'+str(presenceweight)+ 'dw'+ str(dupweight)+ 'll'+str(loss_lambda)+ 'pl'+ str(presenceweight) +'dl' + str(dupl_lamba)
    #print(parastr)
    startdict={'presence':presenceweight, 'loss':lossweight, 'dup':dupweight}
    lambdadict={'presence':presence_lambda, 'loss':loss_lambda, 'dup':dupl_lamba}



    if compile == True:
        with open_file(config_utils.omadir + 'OmaServer.h5', mode="r") as h5_oma:
            lsh_builder = LSHBuilder(h5_oma, saving_folder= config_utils.datadir , saving_name=db, numperm = 256,
            treeweights= None , taxfilter = taxfilter, taxmask= tax_mask , lambdadict= lambdadict, start= startdict)
            hashes, forest , mat = lsh_builder.run_pipeline()
            #hashes, forest, lshpath =lsh_builder.run_pipeline()
        print( 'done compiling')
    else:
        saving_path = config_utils.datadir + db
        hashes = saving_path + 'hashes.h5'
        forest = saving_path + 'newlshforest.pkl'
        mat = saving_path+ 'hogmat.h5'

    print('query DB and calculate error')
    print('load profiler')
    p = profiler.Profiler(lshforestpath = forest, hashes_h5=hashes, mat_path= None )
    print('done')
    print('loading validation')

    folder = config_utils.datadir + 'GOData/'
    if val is None:
        val = validation_semantic_similarity.Validation_semantic_similarity( folder + 'go-basic.obo' ,
            folder + 'goframe.pkl' , folder + 'oma-go.txt' , config_utils.omadir + 'OmaServer.h5' , folder + 'termcounts.pkl' )
    print( 'done')
    print('testing db')

    if not hoglist:
        #sample random hogs
        hoglist = list(np.random.randint(0, high=610000, size=5000, dtype='l'))
        hoglist = [ hashutils.fam2hogid(s)  for s in hoglist]

    scores = {}
    retq = mp.Queue()
    lock = mp.Lock()
    timelimit = 100

    for i,hog in enumerate(hoglist):
        print(hog)
        res = p.hog_query( hog_id = hog , k = 20)
        res = set([ hashutils.fam2hogid(r) for r in res]+[hog])
        #write loop for sem sim check with timeout here
        processes = {}
        for combo in itertools.combinations(res,2):
            processes[combo] = {'time':time.time() , 'process': mp.Process( target = val.semantic_similarity_score_mp , args = (combo[0],combo[1],retq , lock)  ) }
            processes[combo]['process'].start()
            while len(processes)> mp.cpu_count()/4:
                time.sleep(.01)
                for c in processes:
                    if processes[c]['time']>timelimit or processes[c]['process'].exitcode is not None:
                        #print( c[0] +':' + c[1] + ' done')
                        processes[c]['process'].terminate()
                        del(processes[c])
                        break

        while len(processes)> 0:
            time.sleep(.01)
            for c in processes:
                if processes[c]['time']>timelimit or processes[c]['process'].exitcode is not None:
                    processes[c]['process'].terminate()
                    if rocesses[c]['time']>timelimit:
                        print('timeout')
                    del(processes[c])
                    break
        hogsemsim ={}
        while retq.empty() == False:
            combo,semsim = retq.get()
            print(combo)
            print(semsim)

            hogsemsim[combo]=semsim

        scores.update( { combo: {'query_num':i, 'hog_sem_sim': hogsemsim[combo]
        , 'hog_jaccard_sim' : p.hog_v_hog(combo[0],combo[1])
        }  for combo in itertools.combinations(res,2)  if combo in hogsemsim } )
        print(scores)

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

        #try some points
        bo.probe(
        params={'lossweight':  1,
                'presenceweight': 0,
                'dupweight':0,
                'loss_lambda':0,
                'presence_lambda':0,
                'dupl_lamba':0
                },
        lazy=False,
        )

        bo.probe(
        params={'lossweight':  0,
                'presenceweight': 1,
                'dupweight':0,
                'loss_lambda':0,
                'presence_lambda':0,
                'dupl_lamba':0
                },
        lazy=False,
        )

        bo.probe(
        params={'lossweight':  0,
                'presenceweight': 0,
                'dupweight':1,
                'loss_lambda':0,
                'presence_lambda':0,
                'dupl_lamba':0
                },
        lazy=False,
        )

        #intuitively... seems like a good idea
        bo.probe(
        params={'lossweight':  1,
                'presenceweight': .5,
                'dupweight':0,
                'loss_lambda':0,
                'presence_lambda':0,
                'dupl_lamba':0
                },
        lazy=False,
        )

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
