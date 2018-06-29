# should be useless now

import tables
import numpy as np
import pickle
import h5py
import multiprocessing as mp
import pandas as pd
import time
import itertools

from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import TermCounts
from pyoma.browser import db

import config
import validation.phyloValidationGoTerm as validationGoTerm
from profileGen import get_hash_hog_id

# load all the files (oma database, go dag file, association gaf file)
# set oma database
# contains genes and related GO terms
h5OMA = tables.open_file(config.omadir + 'OmaServer.h5', mode='r')
dbObj = db.Database(h5OMA)
omaIdObj = db.OmaIdMapper(dbObj)
# load gene ontology DAG
# use guaranteed acyclic basic GO
# TODO provide correct path for those
go = obo_parser.GODag('./data/go.obo')
associations = read_gaf('./data/gene_association.tair')

# Get the counts of each GO term.
termcounts = TermCounts(go, associations)
goTermAnalysis = validationGoTerm.SemanticSimilarityAnalysis(go, h5OMA, termcounts)

pool = mp.Pool()


def jaccard_mp(h1, h2):
    """
    computes jaccard score for two hashes
    returns score and index
    """
    i, h1 = h1
    j, h2 = h2

    return(i, j, h1.jaccard(h2))


def result_fam_id(result):
    """
    returns the hog fam id for any key of the lsh
    """
    fam_id = str(result.split('-', 1)[0])
    hog_id = "HOG:" + (7-len(fam_id)) * '0' + fam_id

    return hog_id


events = ['duplication', 'gain', 'loss', 'presence']
# contains 2 hog ids, event, go terms, semantic score and jaccard score
result_dict = {}
# only contains the matrices
matrices_dict = {}

with h5py.File(config.datadirLaurent + 'May_09_2018_13_42hashes.h5', 'r') as h5hashes:
    print('loading lsh')
    with open(config.datadirLaurent + 'May_09_2018_13_42_0.7_newlsh.pkl', 'rb') as lsh_file:

        lsh_unpickled = pickle.Unpickler(lsh_file)
        lsh = lsh_unpickled.load()
        print('all loaded')
        # generate random queries from oma
        # gen_rand_queries = 10000
        # np.random.seed(1)
        # queries = list(np.random.randint(low=1, high=len(h5OMA.root.OrthoXML.Index), size=gen_rand_queries))
        # filter the queries; they should contain at least one go terms
        drosoHogList = []
        for hog in dbObj.get_hdf5_handle().get_node('/HogLevel').where('Level == b"Drosophila melanogaster"'):
            drosoHogList.append((hog[0], hog[1]))

        queries =[]

        for fam, hog in drosoHogList:
            if len(queries) > 100:
                break
            # be less restrictive here. Like, only one go term per hog is enough. Create a list with ALL of them.
            # Then pick a few (10, 20, 50 ?) randomly (--> store that in "queries")
            if len(goTermAnalysis._get_go_terms(hog)) > 15:
                queries.append(fam)


        #print('the queries are')
        #print(hogsList)
        #queries = [str(query) for query in hogsList if goTermAnalysis.get_go_terms(query))]
        print(queries)
        # create all combo events, necessary to generate all the hashes
        
        # avoid the double for loop here. Better create a variable containing all the element of the list
        # events_combo = [events_combo for events_combo in itertools.combinations(events, n+1) for n in range(1, len(envents)]
        # maybe there is a function (like itertools.combinations) that does that in a simpler way
        for n in range(1, len(events)):
            for events_combo in itertools.combinations(events, n+1):
                for fam_query in queries:

                    # instead of looping on the double for loop above, maybe better to send all to the function ? 
                    # create a dictionary with key:events_combo value:hash (replace the "get_hash_hog_id" function
                    # instead of fetch (which is buggy) recreate the hashes
                    
                    # get hash for this hog
                    hash_query = get_hash_hog_id(fam_query, h5hashes, events_combo)
                    # get the results for this query in the lsh
                    lsh_results_unfiltered = lsh.query(hash_query)
                    # add the query to the list of results and filter the results
                    # (if no go term, remove hog id)
                    lsh_results_filtered = [fam_query] + [result for result in lsh_results_unfiltered if goTermAnalysis._get_go_terms(result_fam_id(result))]

                    # get a list of hashes, one for each filtered results
                    hashes = [(hog_id, get_hash_hog_id(hog_id, events_combo)) for hog_id in lsh_results_filtered]

                    # prepare the matrix for the semantic distance (go score)
                    semantic_distance = np.zeros((len(lsh_results_filtered), len(lsh_results_filtered)))
                    # double for loop necessary, mp cant be done (yet?)
                    for i, hog1 in enumerate(lsh_results_filtered):
                        for j, hog2 in enumerate(lsh_results_filtered):
                            # get semantic dist between two hog (ids) and save it in matrix and in big dict
                            semantic_dist = goTermAnalysis.semantic_similarity_score(hog1, hog2)
                            semantic_distance[i, j] = semantic_dist
                            result_dict[(hog1, hog2, events_combo)]['semantic_distance'] = semantic_dist
                            # also save the go terms in big dict
                            result_dict[(hog1, hog2, events_combo)]['go terms'] = goTermAnalysis._get_go_terms(hog2)

                    # prepare matrix for the jaccard distances
                    jaccard_distance = np.zeros((len(lsh_results_filtered), len(lsh_results_filtered)))
                    # prepare data for mp stuff
                    hashcompare = [(h1, h2) for h1, h2 in itertools.combinations(enumerate(hashes), 2)]
                    # compute jaccard
                    jaccard_results = pool.map_async(hashcompare, jaccard_mp).get()
                    # feed matrix and big dict with results
                    for jac_res in jaccard_results:
                        i, j, jac_dist = jac_res
                        jaccard_distance[i, j] = jac_dist
                        result_dict[(lsh_results_filtered[i], lsh_results_filtered[j], events_combo)]['jaccard_distance'] = jac_dist

                    # save matrices
                    matrices_dict[(fam_query, events_combo)] = {
                        'semantic_distance': semantic_distance,
                        'jaccard_distance': jaccard_distance}

# transform results in dataframe
results_df = pd.DataFrame.from_dict(result_dict)
# save datafram in a csv file with timestamp
csv_path_filename_date = config.datadirLaurent + 'results' + time.strftime("%Y%m%d-%H%M%S")
results_df.to_csv(csv_path_filename_date, sep='\t')
