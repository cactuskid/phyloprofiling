import _pickle as pickle
import pandas as pd
import h5py
import itertools
import ujson as json
import random

from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import TermCounts

from validation import validation_semantic_similarity
from utils import hashutils

from time import time


class Profiler:

    def __init__(self, lsh_path, hashes_path, obo_file_path, gaf_file_path, h5_go_terms_parents_path):

        self.go_terms_hdf5 = h5py.File(h5_go_terms_parents_path, 'r')
        self.hogs2goterms = self.go_terms_hdf5['hog2goterms']

        self.go = obo_parser.GODag(obo_file_path)
        self.associations = read_gaf(gaf_file_path)

        self.term_counts = TermCounts(self.go, self.associations)
        self.goTermAnalysis = validation_semantic_similarity.Validation_semantic_similarity(self.go,
                                                                                            self.term_counts,
                                                                                            self.go_terms_hdf5)

        lsh_file = open(lsh_path, 'rb')
        lsh_unpickled = pickle.Unpickler(lsh_file)
        self.lsh = lsh_unpickled.load()
        self.hashes = h5py.File(hashes_path, mode='r')

    def hog_query(self, hog_id=None, fam_id=None, events=['duplication', 'gain', 'loss', 'presence'], combination=True):
        """
        Given a hog_id or a fam_id as a query, returns a dictionary containing the results of the LSH.
        :param hog_id: query hog id
        :param fam_id: query fam id
        :param events: list of events one wants to query
        :param combination: Boolean, combination of events or not
        :return: dictionary containing the results of the LSH for the given query
        """
        if hog_id is None and fam_id is None:
            return

        if hog_id is not None:
            fam_id = hashutils.hogid2fam(hog_id)

        # query_hashes dict keys:lminhashname, values:hashes
        # get it from h5hashes instead of recomputing it
        query_dict = {}
        for event in events:
            query_hashe = hashutils.fam2hash_hdf5(fam_id, self.hashes, [event])
            name = str(fam_id) + '-' + event
            query_dict[name] = self.lsh.query(query_hashe)

        if combination:
            for j in range(1, len(events)):
                for i in itertools.combinations(events, j + 1):
                    comb_name = str(fam_id)
                    for array in i:
                        comb_name += '-' + array
                    query_hashe = hashutils.fam2hash_hdf5(fam_id, self.hashes, i)
                    query_dict[comb_name] = self.lsh.query(query_hashe)

        return query_dict

    def results(self, hog_id, fam_id, events, combination, all_vs_all):
        """
        Get jaccard score and semantic value between given hog/id and its results in the LSH
        :param hog_id: hog id
        :param fam_id: fam id
        :param events: evolutionary events
        :param combination: Boolean
        :param all_vs_all: Boolean, True: all vs all, False: query vs results+query
        :return: list of dataframes containing hog ids, semantic value, and jaccard score
        """
        results = self.hog_query(hog_id=hog_id, fam_id=fam_id, events=events, combination=combination)
        see_results(results)

        dataframe_list = []

        for query, results_list in results.items():

            print('getting results for {}'.format(query))
            if all_vs_all:
                results_query_dict = self.results_all_vs_all(query, results_list)
            else:
                results_query_dict = self.results_query(query, results_list)

            results_query_df = pd.DataFrame.from_dict(results_query_dict, orient='index')
            results_query_df['event'] = query

            dataframe_list.append(results_query_df)

            print('{} done'.format(query))

        return dataframe_list

    def filter_results(self, results):
        """
        filters results, remove hogs without GO terms, used for validation
        :param results: dictionary of results
        :return: filtered dictionary of results
        """
        filtered = {}
        # print(results)
        for query, results_list in results.items():

            filtered_list_results = []
            for result in results_list:

                result_fam = hashutils.result2fam(result)
                if result_fam not in filtered_list_results and self.hogs2goterms[result_fam] \
                    and json.loads(self.hogs2goterms[result_fam]):
                    filtered_list_results.append(result)
            filtered[query] = filtered_list_results
        return filtered

    def results_all_vs_all(self, query, results_list):
        results_dict = {}
        results_list = [query] + results_list
        for hog_event_1 in results_list:

            for hog_event_2 in results_list:
                results_dict.update(self.get_scores(hog_event_1, hog_event_2, results_dict))
        return results_dict

    def results_query(self, query, results_list):
        results_dict = {}
        hog_event_1 = query
        results_list = [query] + results_list

        for hog_event_2 in results_list:
                results_dict.update(self.get_scores(hog_event_1, hog_event_2, results_dict))

        return results_dict

    def get_scores(self, hog_event_1, hog_event_2, results_dict):
        hog1 = hashutils.result2hogid(hog_event_1)
        hog2 = hashutils.result2hogid(hog_event_2)

        semantic_dist = self.compute_semantic_distance(hog1, hog2)
        results_dict[(hog1, hog2)] = {'Semantic': semantic_dist}
        jaccard_dist = self.compute_jaccard_distance(hog_event_1, hog_event_2)
        results_dict[(hog1, hog2)]['Jaccard'] = jaccard_dist

        return results_dict


    # TODO: add get string functions hashes_error_files

    def compute_string_score(self):
        pass

    def hog2string(self):
        pass

    def string2interactions(self):
        pass

    def compute_semantic_distance(self, hog_1, hog_2):

        semantic_dist = self.goTermAnalysis.semantic_similarity_score(hog_1, hog_2)

        return semantic_dist

    def compute_jaccard_distance(self, hog_event_1, hog_event_2):

        hash_1, hash_2 = [hashutils.fam2hash_hdf5(hashutils.result2fam(hog_event), self.hashes)
                          for hog_event in [hog_event_1, hog_event_2]]
        jaccard_dist = hash_1.jaccard(hash_2)

        return jaccard_dist

    def get_hogs_with_annotations(self):
        print('getting hogs')
        hogs_with_annotations = []
        for fam, goterms in enumerate(self.hogs2goterms):
            if goterms:
                if json.loads(goterms):
                    hogs_with_annotations.append(fam)

        print(len(hogs_with_annotations))
        return hogs_with_annotations

    def validate_pipeline_go_terms(self, path_to_save):

        hogs_with_annotations = self.get_hogs_with_annotations()
        # for each hog with annotations, query results

        dataframe_list = []

        for hog in hogs_with_annotations:
            raw_results = self.hog_query(fam_id=hog)
            filtered_results = self.filter_results(raw_results)
            see_results(filtered_results)
            for query, results_list in filtered_results.items():
                start_time = time()
                results_query_dict = self.results_query(query, results_list)
                results_query_df = pd.DataFrame.from_dict(results_query_dict, orient='index')
                results_query_df['event'] = query

                dataframe_list.append(results_query_df)

                print('{} done in {}'.format(query, time() - start_time))

        concat_result = pd.concat(list(dataframe_list))
        concat_result.to_csv(path_to_save, sep='\t')
        print('DONE')

    def validate_pipeline(self, path_to_hog_id_file, path_to_save):

        # get list of queries
        hog_ids = list(pd.read_csv(path_to_hog_id_file, sep='\t')['hog id'])

        dataframe_list = []

        # query the LSH, get dict query:list of result
        for hog in hog_ids:
            raw_results = self.hog_query(fam_id=hog)
            filtered_results = self.filter_results(raw_results)
            see_results(filtered_results)

            for query, list_results in filtered_results.items():
                start_time = time()

                number_samples = 10 if len(list_results) >= 10 else len(list_results)
                random_results = random.sample(list_results, number_samples)

                results_query_dict = self.results_all_vs_all(query, random_results)

                results_query_df = pd.DataFrame.from_dict(results_query_dict, orient='index')
                results_query_df['event'] = query

                dataframe_list.append(results_query_df)

                print('{} done in {}'.format(query, time() - start_time))

        concat_result = pd.concat(list(dataframe_list))
        concat_result.to_csv(path_to_save, sep='\t')
        print('DONE')

    def validate_pipeline_string(self, path_to_save):



    def query_pipeline(self, hog_id=None, fam_id=None, events=['duplication', 'gain', 'loss', 'presence'],
                     combination=True, path_to_save=None, all_vs_all=False):

        print('getting results')
        df_results = self.results(hog_id=hog_id, fam_id=fam_id, events=events, combination=combination,
                                  all_vs_all=all_vs_all)
        print('saving ...')
        if path_to_save is not None:
            concat_result = pd.concat(list(df_results))
            concat_result.to_csv(path_to_save, sep='\t')
            print(df_results)
        print('done')


def see_results(results):
    for k, v in results.items():
        print('{} queries found for {}'.format(len(v), k))
