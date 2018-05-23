import _pickle as pickle
import tables
import pandas as pd
import multiprocessing as mp
import itertools
import h5py

from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import TermCounts
from pyoma.browser import db

from validation import phyloValidationGoTerm as validationGoTerm
from utils import files_utils, hashutils, config_utils

# to remove
from time import time


class Profiler:

    def __init__(self, h5_oma_path):
        # load all the files (oma database, go dag file, association gaf file=
        self.h5OMA = tables.open_file(h5_oma_path, mode='r')
        self.db_obj = db.Database(self.h5OMA)
        self.omaIdObj = db.OmaIdMapper(self.db_obj)
        self.replacement_dic, self.tree = files_utils.get_species_tree_replacement_dic(self.h5OMA, self.omaIdObj)

        self.go = None
        self.associations = None
        self.term_counts = None
        self.goTermAnalysis = None

        self.lsh = None
        # init here, so if multiple query with similar hog id, their go terms are already loaded and computed

        self.best_go_terms_hog_id_dict = {}

    def go_benchmarking_init(self, obo_file_path, gaf_file_path):
        self.go = obo_parser.GODag(obo_file_path)
        self.associations = read_gaf(gaf_file_path)
        # Get the counts of each GO term.
        self.term_counts = TermCounts(self.go, self.associations)
        self.goTermAnalysis = validationGoTerm.SemanticSimilarityAnalysis(self.go, self.h5OMA, self.term_counts)

    def lsh_loader(self, lsh_path):
        lsh_file = open(lsh_path, 'rb')
        lsh_unpickled = pickle.Unpickler(lsh_file)
        self.lsh = lsh_unpickled.load()

    def hog_query(self, hog_id=None, fam_id=None, events=['duplication', 'gain', 'loss', 'presence'], combination=True):

        if hog_id is None and fam_id is None:
            return

        if hog_id is not None:
            fam_id = hashutils.hogid2fam(hog_id)

        # query_hashes dict keys:lminhashname, values:hashes
        query_hashes = hashutils.fam2hashes(fam_id, self.db_obj, self.tree, self.replacement_dic, events, combination)

        query_dict = {}
        # query_hashes['dict'] contains all 15 leanminhashes
        for name, hashvalue in query_hashes['dict'].items():
            # lsh.query returns a list of hash (fam-event1-event2- etc)
            query_dict[name] = self.lsh.query(hashvalue)

        return query_dict

    def results(self, hog_id=None, fam_id=None, events=['duplication', 'gain', 'loss', 'presence'],
                combination=True, scores=False):

        print('getting queries')
        results = self.hog_query(hog_id=hog_id, fam_id=fam_id, events=events, combination=combination)
        for k, v in results.items():
            print('{} queries found for {}'.format(len(v), k))
        if not scores:
            results_df = pd.DataFrame.from_dict(results)
            return results_df

        start_time = time()
        self.get_best_go_term_per_hog_id(results)
        print('time to load best got : {}'.format(time()-start_time))

        result_dict = {}
        start_time = time()
        # semantic scores
        for query, result in results.items():
            if 'duplication-loss-presence-gain' in query:
                print('computing scores for query {}'.format({query}))

                related_hog_list = list(set([query] + [r for r in result]))
                events_combo = hashutils.result2events(query)
                # result_dict = self.compute_semantic_distance(related_hog_list, query, result_dict, events_combo)

                # get jaccard score, here its possible to use mp
                result_dict = self.compute_jaccard_distance(related_hog_list, query, result_dict, events_combo)
                print('computing scores for query {} took {}'.format(query, time()-start_time))

        print(result_dict)
        results_df = pd.DataFrame.from_dict(result_dict, orient='index')
        return results_df

    def get_best_go_term_per_hog_id(self, results):
        # new function !
        # best_go_term dictonary updater

        hogs_id_list = []

        for query, result in results.items():
            result_list = list(set([query] + [r for r in result]))
            result_hogs_id_list = list(set([hashutils.result2hogid(result) for result in result_list]))
            hogs_id_list = list(set(hogs_id_list + result_hogs_id_list))

        start_time_big = time()
        for hog in hogs_id_list:
            start_time = time()
            if hog not in self.best_go_terms_hog_id_dict.keys() and len(self.best_go_terms_hog_id_dict.keys())<20:
                print('looking for go terms for {}'.format(hog))
                self.best_go_terms_hog_id_dict[hog] = self.goTermAnalysis.get_best_go_terms_per_hog_id(hog)
                print(time()-start_time)

        print('time to update the dict {}'.format(time()-start_time_big))

    def compute_semantic_distance(self, hogs_list, query, result_dict, events_combo):

        events_combo_string = ''.join(event for event in sorted(events_combo))

        for i, hog_event_1 in enumerate(hogs_list):
            start_time = time()

            for j, hog_event_2 in enumerate(hogs_list):
                hog1 = hashutils.result2hogid(hog_event_1)
                hog2 = hashutils.result2hogid(hog_event_2)
                if hog1 in self.best_go_terms_hog_id_dict.keys() and hog2 in self.best_go_terms_hog_id_dict.keys():
                    start_time = time()

                    # get semantic dist between two hog (ids) and save it in matrix and in big dict
                    semantic_dist = self.goTermAnalysis.semantic_similarity_score_from_go_terms(self.best_go_terms_hog_id_dict[hog1], self.best_go_terms_hog_id_dict[hog2])
                    result_dict[(hog1, hog2, events_combo_string, 'semantic_distance')] = semantic_dist
                    # also save the go terms in big dict

                    # result_dict[(query, hog1, hog2, events_combo_string, 'go term hog 1')] = self.goTermAnalysis.get_go_terms(hog1)
                    # result_dict[(query, hog1, hog2, events_combo_string, 'go term hog 2')] = self.goTermAnalysis.get_go_terms(hog2)
            print('computing results for {}, took {}'.format(hog_event_1, time()-start_time))

        return result_dict

    def compute_jaccard_distance(self, hogs_list, query, result_dict, events_combo):
        # prepare data for mp stuff
        print('start computing jaccard')


        events_combo_string = ''.join(event for event in sorted(events_combo))
        print('before hashes')
        start = time()

        with h5py.File(config_utils.datadirLaurent + 'May_16_2018_16_07hashes.h5', 'r') as h5hashes:

            hashes = {result: hashutils.get_hash_hog_id(fam=hashutils.result2fam(result), h5mat=h5hashes, events=events_combo) for result in hogs_list}

        print(list(hashes.values())[0])
        # take from hash5 instead
        print('after hashes {}'.format(time()-start))
        hash_compare = [((h1[0],hashes[h1[1]]), (h2[0],hashes[h2[1]])) for h1, h2 in itertools.combinations(enumerate(hashes), 2)]
        # compute jaccard
        pool = mp.Pool(10)
        print('before jac')
        jaccard_results = pool.map_async(jaccard_mp, hash_compare, 20).get()


        print(jaccard_results)
        # feed matrix and big dict with results
        print('after jac')
        for jac_res in jaccard_results:

            i, j, jac_dist = jac_res
            hog1 = hashutils.result2hogid(hogs_list[i])
            hog2 = hashutils.result2hogid(hogs_list[j])

            result_dict[(hog1, hog2, events_combo_string, 'jaccard_distance')] = jac_dist #
            print(result_dict)
        pool.close()

        return result_dict

    def results_save(self, hog_id=None, fam_id=None, events=['duplication', 'gain', 'loss', 'presence'],
                     combination=True, scores=False, path_to_save=None):
        print('getting results')
        df_results = self.results(hog_id=hog_id, fam_id=fam_id, events=events, combination=combination, scores=scores)
        print('saving ...')
        if path_to_save is not None:
            df_results.to_csv(path_to_save, sep='\t')
            print(df_results)
        print('done')


def jaccard_mp(args):
    """
    computes jaccard score for two hashes
    returns score and index
    """
    hash1, hash2 = args
    i, h1 = hash1
    j, h2 = hash2

    print(type(h1))
    print(h1)

    return i, j, h1.jaccard(h2)
