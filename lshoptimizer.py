import lsh_builder
import profiler
import validation_semantic_similarity

import config_utils
import pickle

from bayes_opt import BayesianOptimization


def profiling_error( db , taxfilter, tax_mask lossweight , presenceweight, dupweight, loss_lambda , presence_lambda , dupl_lamba, presenceweight ):


    #compile lsh
    with open_file(config_utils.omadir + 'OmaServer.h5', mode="r") as h5_oma:
        lsh_builder = LSHBuilder(h5_oma, saving_folder= config_utils.datadir , saving_name=db, numperm = 256,
        treeweights= None , taxfilter = taxfilter, taxmask= tax_mask , lambdadict= lambdadict, start= startdict)
        lsh, forest, hashes = lsh_builder.run_pipeline()
    #init profiler
    p = profiler.Profiler(lsh_path = forest, hashes_path = hashes, oma_path= config_utils.omadir + 'OmaServer.h5', mat_path = None,
    unimap_path = None, string_data_path = None, taxfilter = None, taxmast = None, weights = None, GO= )

    #
    if randomHOGS = False:
        #gold standard list
        rando = list(np.randint())
        results = [ p.hog_query(r) for r in rando ]
        godists = [ ]
        go_enrich = [ ]
    else:
        # 1000 random Hogs with annotation
        p()

    return 1/errorval


if __name__ == 'main':
    dbdict = {
    'plants': { 'taxfilter': None , 'taxmask': 33090 }
    'all': { 'taxfilter': None , 'taxmask': None }
    'archaea':{ 'taxfilter': None , 'taxmask': 2157 }
    'bacteria':{ 'taxfilter': None , 'taxmask': 2 }
    'eukarya':{ 'taxfilter': None , 'taxmask': 2759 }
    'protists':{ 'taxfilter': [2 , 2157 , 33090 , 4751, 33208] , 'taxmask':None }
    'fungi':{ 'taxfilter': None , 'taxmask': 4751 }
    'metazoa':{ 'taxfilter': None , 'taxmask': 33208 }
    }
    for db in dbdict:

        error = partial()

        #get error for the first point with all weights at 1
        
        bo = BayesianOptimization( profiling_error ,  {'lossweight': (0, 1),
                                                        'presencweight': (0, 1),
                                                        'dupweight':(0,1),
                                                        'loss_lambda':(0,1),
                                                        'presence_lambda':(0,1),
                                                        'dupl_lamba':(0,1)
                                                        })

        # Hit the big red button.
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
