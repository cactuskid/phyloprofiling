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
from utils import files_utils, hashutils


class Profiler:

    def __init__(self, h5_oma_path):
        self.h5OMA = tables.open_file(h5_oma_path, mode='r')
        self.db_obj = db.Database(self.h5OMA)
        self.omaIdObj = db.OmaIdMapper(self.db_obj)
        self.replacement_dic, self.tree = files_utils.get_species_tree_replacement_dic(self.h5OMA, self.omaIdObj)

        self.go = None
        self.associations = None
        self.term_counts = None
        self.goTermAnalysis = None

        self.lsh = None
        self.hashes = None

        self.go_terms_dict = {}

    def go_benchmarking_init(self, obo_file_path, gaf_file_path):
        self.go = obo_parser.GODag(obo_file_path)
        self.associations = read_gaf(gaf_file_path)
        # Get the counts of each GO term.
        self.term_counts = TermCounts(self.go, self.associations)
        self.goTermAnalysis = validationGoTerm.SemanticSimilarityAnalysis(self.go, self.h5OMA, self.term_counts)

    def lsh_loader(self, lsh_path, hashes_path):
        lsh_file = open(lsh_path, 'rb')
        lsh_unpickled = pickle.Unpickler(lsh_file)
        self.lsh = lsh_unpickled.load()

        self.hashes = hashes_path

    def hog_query(self, hog_id=None, fam_id=None, events=['duplication', 'gain', 'loss', 'presence'], combination=True):

        if hog_id is None and fam_id is None:
            return

        if hog_id is not None:
            fam_id = hashutils.hogid2fam(hog_id)

        # query_hashes dict keys:lminhashname, values:hashes
        # get it from h5hashes instead of recomputing it
        query_hashes = hashutils.fam2hashes(fam_id, self.db_obj, self.tree, self.replacement_dic, events, combination)

        query_dict = {}
        # query_hashes['dict'] contains all 15 leanminhashes
        for name, hashvalue in query_hashes['dict'].items():
            # lsh.query returns a list of hash (fam-event1-event2- etc)
            query_dict[name] = self.lsh.query(hashvalue)

        return query_dict

    def validation(self, hog_id=None, fam_id=None, events=['duplication', 'gain', 'loss', 'presence'],
                   combination=True, scores=True):

        results = self.hog_query(hog_id=hog_id, fam_id=fam_id, events=events, combination=combination)
        for k, v in results.items():
            print('{} queries found for {}'.format(len(v), k))
            print(v)
        if not scores:
            results_df = pd.DataFrame.from_dict(results)
            return results_df

        validation_list = []

        for query, results_list in results.items():
            if 'duplication' not in query and 'presence' in query:
                validation_dict = {}
                print('getting score for {}'.format(query))

                validation_dict.update(self.compute_jaccard_distance_query(query, results_list, validation_dict, events))

                validation_dict.update(self.compute_semantic_distance_query(results_list, validation_dict, events, query))

                validation_df = pd.DataFrame.from_dict(validation_dict, orient='index')
                validation_df['event'] = query

                print(validation_df)

                validation_list.append(validation_df)

        return validation_list

    def results(self, hog_id=None, fam_id=None, events=['duplication', 'gain', 'loss', 'presence'],
                combination=True, scores=False):
        results = self.hog_query(hog_id=hog_id, fam_id=fam_id, events=events, combination=combination)
        for k, v in results.items():
            print('{} queries found for {}'.format(len(v), k))
            print(v)
        if not scores:
            results_df = pd.DataFrame.from_dict(results)
            return results_df

        dataframe_list = []

        # semantic scores
        for query, result in results.items():
            result_dict = {}
            print('computing scores for query {}'.format({query}))

            related_hog_list = list(set([query] + [r for r in result]))
            events_combo = hashutils.result2events(query)

            result_dict.update(self.compute_semantic_distance(related_hog_list, result_dict, events_combo))

            # get jaccard score, here its possible to use mp
            jaccard = self.compute_jaccard_distance(related_hog_list, result_dict, events_combo)
            for index in jaccard:
                result_dict[index].update(jaccard[index])

            results_df = pd.DataFrame.from_dict(result_dict, orient='index')
            dataframe_list.append(results_df)

        return dataframe_list

    def compute_semantic_distance(self, hogs_list, result_dict, events_combo):

        events_combo_string = ''.join(event for event in sorted(events_combo))
        for i, hog_event_1 in enumerate(hogs_list):
            for j, hog_event_2 in enumerate(hogs_list):
                hog1 = hashutils.result2hogid(hog_event_1)
                hog2 = hashutils.result2hogid(hog_event_2)


                semantic_dist = self.goTermAnalysis.semantic_similarity_score(hog1, hog2)
                result_dict[(hog1, hog2)] = {'Semantic': semantic_dist}

        return result_dict

    def compute_semantic_distance_query(self, hogs_list, result_dict, events_combo, query):

        hog_query = hashutils.result2hogid(query)

        if hog_query not in self.go_terms_dict:
            self.go_terms_dict[hog_query] = self.goTermAnalysis.get_go_terms(hog_query)
            print(self.go_terms_dict)

        for i, hog_event in enumerate(hogs_list):
            hog_other = hashutils.result2hogid(hog_event)

            if hog_other not in self.go_terms_dict:
                self.go_terms_dict[hog_other] = self.goTermAnalysis.get_go_terms(hog_other)

            semantic_dist = self.goTermAnalysis.semantic_similarity_score_from_go_terms(self.go_terms_dict[hog_query],
                                                                                        self.go_terms_dict[hog_other])
            result_dict[(hog_query, hog_other)]['Semantic'] = semantic_dist

        return result_dict

    def compute_jaccard_distance_query(self, query, results_list, result_dict, events_combo):

        events_combo_string = ''.join(event for event in sorted(events_combo))

        with h5py.File(self.hashes, 'r') as h5hashes:
            query_hash = {query: hashutils.fam2hash_hdf5(hashutils.result2fam(query), h5hashes, events_combo)}

            results_hashes = {result: hashutils.fam2hash_hdf5(hashutils.result2fam(result), h5hashes, events_combo)
                              for result in results_list}

        print(results_hashes)

        hash_compare_query = [((0, query_hash[query]), (i, results_hashes[h])) for i, h in enumerate(results_hashes)]
        print(hash_compare_query)
        #pool = mp.Pool(10)
        print('before pool')
        #jaccard_results = pool.map_async(jaccard_mp, hash_compare_query, 10).get(9999999)

        jaccard_results = list(map(jaccard_mp, hash_compare_query))

        # feed matrix and big dict with results
        for jac_res in jaccard_results:
            print('getting results')
            i, j, jac_dist = jac_res
            hog1 = hashutils.result2hogid(query)
            hog2 = hashutils.result2hogid(results_list[j])

            result_dict[(hog1, hog2)] = {'Jaccard': jac_dist}

        #pool.close()

        return result_dict

    def compute_jaccard_distance(self, hogs_list, result_dict, events_combo):

        events_combo_string = ''.join(event for event in sorted(events_combo))

        with h5py.File(self.hashes, 'r') as h5hashes:

            hashes = {result: hashutils.fam2hash_hdf5(hashutils.result2fam(result), h5hashes, events_combo)
                      for result in hogs_list}


        hash_compare = [((h1[0], hashes[h1[1]]), (h2[0], hashes[h2[1]]))
                        for h1, h2 in itertools.combinations(enumerate(hashes), 2)]

        pool = mp.Pool(10)

        jaccard_results = pool.map_async(jaccard_mp, hash_compare, 20).get()

        # feed matrix and big dict with results
        for jac_res in jaccard_results:

            i, j, jac_dist = jac_res
            hog1 = hashutils.result2hogid(hogs_list[i])
            hog2 = hashutils.result2hogid(hogs_list[j])

            result_dict[(hog1, hog2)] = {'Jaccard': jac_dist}

        pool.close()

        return result_dict

    def results_save(self, hog_id=None, fam_id=None, events=['duplication', 'gain', 'loss', 'presence'],
                     combination=True, scores=False, path_to_save=None):
        print('getting results')
        df_results = self.results(hog_id=hog_id, fam_id=fam_id, events=events, combination=combination, scores=scores)
        print('saving ...')
        if path_to_save is not None:

            concat_result = pd.concat(list(df_results))
            concat_result.to_csv(path_to_save, sep='\t')
            print(df_results)
        print('done')

    def validation_save(self, hog_id=None, fam_id=None, events=['duplication', 'gain', 'loss', 'presence'],
                        combination=True, scores=True, path_to_save=None):

        print('getting results')
        df_results = self.validation(hog_id=hog_id, fam_id=fam_id, events=events, combination=combination,
                                  scores=scores)
        print('saving ...')
        if path_to_save is not None:
            concat_validation = pd.concat(df_results)
            concat_validation.to_csv(path_to_save, sep='\t')
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

    return i, j, h1.jaccard(h2)
