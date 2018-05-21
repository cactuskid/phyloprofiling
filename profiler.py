import pickle
import tables
import pandas as pd
import multiprocessing as mp
import itertools

from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import TermCounts
from pyoma.browser import db

from validation import phyloValidationGoTerm as validationGoTerm
from utils import files_utils, hashutils

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
        for k,v in results.items():
            print('{} queries found for {}'.format(len(v),k))
        if not scores:
            results_df = pd.DataFrame.from_dict(results)
            return results_df



        self.get_best_go_term_per_hog_id(results)

        result_dict = {}

        # semantic scores
        for query, result in results.items():
            print('computing scores for query {}'.format({query}))
            related_hog_list = list(set([query] + [r for r in result]))
            events_combo = hashutils.result2events(query)
            result_dict = self.compute_semantic_distance(related_hog_list, query, result_dict, events_combo)

            # get jaccard score, here its possible to use up
            result_dict = self.compute_jaccard_distance(related_hog_list, query, result_dict, events_combo)

        return result_dict

    def get_best_go_term_per_hog_id(self, results):
        # new function !
        # best_go_term dictonary updater

        hogs_id_list = []

        for query, result in results.items():
            result_list = list(set([query] + [r for r in result]))
            hogs_id_list = list(set([hashutils.result2hogid(hog) for hog in result_list]))

        for hog in hogs_id_list:
            if hog not in self.best_go_terms_hog_id_dict.keys():
                self.best_go_terms_hog_id_dict[hog] = self.goTermAnalysis.API_get_best_go_terms_per_gene(hog)

    def compute_semantic_distance(self, hogs_list, query, result_dict, events_combo):

        events_combo_string = ''.join(event for event in sorted(events_combo))

        for i, hog_event_1 in enumerate(hogs_list):
            for j, hog_event_2 in enumerate(hogs_list):
                hog1 = hashutils.result2hogid(hog_event_1)
                hog2 = hashutils.result2hogid(hog_event_2)

                # get semantic dist between two hog (ids) and save it in matrix and in big dict
                semantic_dist = self.goTermAnalysis.API_compute_score_method_2.(self.best_go_terms_hog_id_dict[hog1], self.best_go_terms_hog_id_dict[hog2])
                result_dict[(query, hog1, hog2, events_combo_string, 'semantic_distance')] = semantic_dist
                # also save the go terms in big dict

                #result_dict[(query, hog1, hog2, events_combo_string, 'go term hog 1')] = self.goTermAnalysis.get_go_terms(hog1)
                #result_dict[(query, hog1, hog2, events_combo_string, 'go term hog 2')] = self.goTermAnalysis.get_go_terms(hog2)


        return result_dict

    def jaccard_mp(self, hash1, hash2):
        """
        computes jaccard score for two hashes
        returns score and index
        """
        i, h1 = hash1
        j, h2 = hash2

        return i, j, h1.jaccard(h2)

    def compute_jaccard_distance(self, hogs_list, query, result_dict, events_combo):
        # prepare data for mp stuff
        pool = mp.Pool()

        events_combo_string = ''.join(event for event in sorted(events_combo))
        hashes = {result: hashutils.result2hash(result, self.tree, self.replacement_dic, self.db_obj) for result in hogs_list}

        hash_compare = [(h1, h2) for h1, h2 in itertools.combinations(enumerate(hashes), 2)]
        # compute jaccard
        jaccard_results = pool.map_async(hash_compare, self.jaccard_mp).get()
        # feed matrix and big dict with results
        for jac_res in jaccard_results:
            i, j, jac_dist = jac_res
            hog1 = hashutils.result2hogid(hogs_list[i])
            hog2 = hashutils.result2hogid(hogs_list[j])
            result_dict[(query, hog1, hog2, events_combo_string, 'jaccard_distance')] = jac_dist

        return result_dict

    def results_save(self, hog_id=None, fam_id=None, events=['duplication', 'gain', 'loss', 'presence'],
                     combination=True, scores=False, path_to_save=None):
        print('getting results')
        df_results = self.results(hog_id=hog_id, fam_id=fam_id, events=events, combination=combination, scores=scores)
        print('saving ...')
        if path_to_save is not None:
            df_results.to_csv(path_to_save, sep='\t')
        print('done')