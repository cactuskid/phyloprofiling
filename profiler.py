import _pickle as pickle
import pandas as pd
import h5py
import itertools
import ujson as json
import random
from scipy.sparse import csr_matrix

from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import TermCounts
from pyoma.browser import db

from validation import validation_semantic_similarity
from utils import hashutils, string_utils


from time import time
class Profiler:

    def __init__(self, oma_path ,lsh_path, hashes_path, mat_path = None , unimap_path = None ,string_data_path = None , taxfilter = None, taxmast = None, weights = None ):
        #use the lsh forest or the lsh


        self.go_terms_hdf5 = h5py.File(h5_go_terms_parents_path, 'r')
        self.hogs2goterms = self.go_terms_hdf5['hog2goterms']

        self.go = obo_parser.GODag(obo_file_path)
        self.associations = read_gaf(gaf_file_path)

        self.term_counts = TermCounts(self.go, self.associations)
        self.goTermAnalysis = validation_semantic_similarity.Validation_semantic_similarity(self.go,
                                                                                            self.term_counts,
                                                                                            self.go_terms_hdf5)

        self.h5OMA = oma_path
        self.db_obj = db.Database(self.h5OMA)

        lsh_file = open(lsh_path, 'rb')
        lsh_unpickled = pickle.Unpickler(lsh_file)
        self.lsh = lsh_unpickled.load()
        self.hashes_h5 = h5py.File(hashes_path, mode='r')

        if unimap_path:
            self.unimap_h5 = h5py.File(unimap_h5, mode='r')

        if string_data_path:
            self.string_data_path = string_data_path
            self.r1 = string_stringdataMap.connect2IDmap()
            self.r2 = string_stringdataMap.connect2Stringmap()

        if mat_path:
            profile_matrix_file = open(profile_matrix_path, 'rb')
            profile_matrix_unpickled = pickle.Unpickler(profile_matrix_file)
            self.profile_matrix = profile_matrix_unpickled.load()


    def hog_query(self, hog_id=None, fam_id=None):
        """
        Given a hog_id or a fam_id as a query, returns a dictionary containing the results of the LSH.
        :param hog_id: query hog id
        :param fam_id: query fam id
        :return: list containing the results of the LSH for the given query
        """
        if hog_id is None and fam_id is None:
            return
        if hog_id is not None:
            fam_id = hashutils.hogid2fam(hog_id)
        # query_hashes dict keys:lminhashname, values:hashes
        # get it from h5hashes instead of recomputing it
        query_dict = {}
        query_hash = hashutils.fam2hash_hdf5_multiset(fam_id, self.hashes_h5)
        results = self.lsh.query(query_hashe)
        return query_dict

    def pull_hashes(self , hoglist):
        return { hog:hashutils.fam2hash_hdf5(hog, self.hashes_h5 ) for hog in hoglist}

    """
    def pull_mapping(self, hoglist):
        #grab the crossrefs for a list of hogs
        return { hog:{ dataset: json.loads(unimap_h5[dataset][hog]) for dataset in unimap_h5 }  for hog in hoglist }
"""
    def pull_go(self,hoglist):
        pass


    def pull_matrows(fams):
        """
        given a list of fams return the submatrix containing their profiles

        :return:fams sorted, sparse mat
        """
        return self.profile_matrix[fams,:]

    def sort_hashes(query_hash,hashes):
        jaccard=[ query_hash.jaccard(hashes[hog]) for hog in hashes]
        index = np.argsort(jaccard)
        sortedhogs = np.asarry(list(hashes.keys()))[index]
        jaccard= jaccard[index]
        return sortedhogs, jaccard

    def allvall_hashes(hashes):
        hashmat = np.zeros((len(hashes),len(hashes)))
        for i , hog1 in enumerate(hashes):
            for j, hog2 in enumerate(hahes):
                hashmat[i,j]= hashes[hog1].jaccard(hashes[hog2])
        return hashmat

    def get_vpairs(fam):
        #get pairwise distance matrix of OMA all v all
        ## TODO: unfinished

        taxa = self.db_obj.hog_levels_of_fam(fam)
        subtaxindex =  { tup[1]:tup[0] for tup in enumerate(taxa)}
        prots = self.db_obj.hog_members_from_hog_id(fam,  'LUCA')
        for prot in prots:
            taxon = prot.ncbi_taxon_id()
            pairs = self.db_obj.get_vpairs(prot)
            for EntryNr1, EntryNr2, RelType , score , distance in list(pairs):
                pass

        return sparsemat , densemat
"""

class Validation_Profiler:

    def __init__(self, lsh_path, hashes_path, obo_file_path, gaf_file_path, h5_go_terms_parents_path, oma_path, string_data_path=None, mat_path=None):
        Profiler.__init__(lsh_path, hashes_path, mat_path)
        self.go_terms_hdf5 = h5py.File(h5_go_terms_parents_path, 'r')
        self.hogs2goterms = self.go_terms_hdf5['hog2goterms']

        self.go = obo_parser.GODag(obo_file_path)
        self.associations = read_gaf(gaf_file_path)

        self.term_counts = TermCounts(self.go, self.associations)
        self.goTermAnalysis = validation_semantic_similarity.Validation_semantic_similarity(self.go,
                                                                                            self.term_counts,
                                                                                            self.go_terms_hdf5)


        self.h5OMA = oma_path
        self.db_obj = db.Database(self.h5OMA)


    def results_query(self, query, results_list):
        results_dict = {}
        hog_event_1 = query
        results_list = [query] + results_list
        for hog_event_2 in results_list:
                results_dict.update(self.get_scores(hog_event_1, hog_event_2, results_dict))
        return results_dict



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


    def get_scores(self, hog_event_1, hog_event_2, results_dict):
        hog1 = hashutils.result2hogid(hog_event_1)
        hog2 = hashutils.result2hogid(hog_event_2)

        semantic_dist = self.compute_semantic_distance(hog1, hog2)
        results_dict[(hog1, hog2)] = {'Semantic': semantic_dist}
        jaccard_dist = self.compute_jaccard_distance(hog_event_1, hog_event_2)
        results_dict[(hog1, hog2)]['Jaccard'] = jaccard_dist

        return results_dict

    def results_with_string_score(self, query, results_list):
        results_dict = {}
        results_list = [query] + results_list
        for hog_event_1 in results_list:
            for hog_event_2 in results_list:

                if len(results_dict) < 10:
                    results = self.get_string_jaccard_scores(hog_event_1, hog_event_2)
                    if results:
                        results_dict.update(results)
        return results_dict

    def get_string_jaccard_scores(self, hog_event_1, hog_event_2):
        hog1 = hashutils.result2hogid(hog_event_1)
        hog2 = hashutils.result2hogid(hog_event_2)

        # TODO get string results list
        string_results_list = self.get_string_scores(hog1, hog2)

        results_dict = {}

        if string_results_list:
            results_dict[(hog1, hog2)] = {'String': string_results_list}

            jaccard_dist = self.compute_jaccard_distance(hog_event_1, hog_event_2)
            results_dict[(hog1, hog2)]['Jaccard'] = jaccard_dist

        return results_dict

    # TODO: add get string functions hashes_error_files

    def get_string_scores(self, hog1, hog2):

        allstring1 = string_stringdataMap.fam2stringID(self.db_obj, hog1, self.r1)
        allstring2 = string_stringdataMap.fam2stringID(self.db_obj, hog2, self.r1)

        if len(allstring1) > 0 and len(allstring2) > 0:

            # string_results = string_stringdataMap.HOGvsHOG(allstring1, allstring2, self.r2, self.string_data_path)
            # here, change function ...
            string_results = string_utils.get_interactions(allstring1, allstring2)
            print(string_results)
            return string_results
        else:
            return None

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

    def validate_pipeline_go_terms(self, path_to_hog_id_file, path_to_save):

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

        # get randomly n hogs
        hog_ids = self.get_random_hogs_with_string_id(2)

        dataframe_list = []

        for hog in hog_ids:
            results = self.hog_query(fam_id=hog)
            see_results(results)

            for query, list_results in results.items():
                start_time = time()

                # do not take all results, otherwise too long
                # number_samples = 10 if len(list_results) >= 10 else len(list_results)
                # random_results = random.sample(list_results, number_samples)

                results_query_dict = self.results_with_string_score(query, list_results)

                results_query_df = pd.DataFrame.from_dict(results_query_dict, orient='index')
                results_query_df['event'] = query

                dataframe_list.append(results_query_df)

                print('{} done in {}'.format(query, time() - start_time))

        concat_result = pd.concat(list(dataframe_list))
        concat_result.to_csv(path_to_save, sep='\t')
        print('DONE')

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

    def query_analysis_pipeline(self, hog_id=None, fam_id=None, events=['duplication', 'gain', 'loss', 'presence'],
                                combination=True, path_to_save=None, all_vs_all=False):

        print('getting results')
        results = self.hog_query(hog_id=hog_id, fam_id=fam_id, events=events, combination=combination)
        see_results(results)

        results_ids = get_hog_ids_from_results(results)

        # get X (matrix)
        results_matrix = self.profile_matrix

    def get_submatrix_form_results(self, results):

        res_mat_list = []

        for query, result in results.items():
            res_mat = csr_matrix((len(result), self.profile_matrix.shape[1]))
            for i, r in enumerate(result):
                res_mat[i, :] = self.profile_matrix[r, :]

            res_mat_list.append(res_mat)

        return res_mat_list

    def get_random_hogs_with_string_id(self, number_hogs):
        # TODO put correct number of hogs in OMA
        hogs_in_OMA = 500000

        fam_ids = []
        while len(fam_ids) < number_hogs:
            rand_hog = random.randint(1, hogs_in_OMA)
            if rand_hog not in fam_ids:
                string = string_stringdataMap.fam2stringID(self.db_obj, hashutils.fam2hogid(rand_hog), self.r1)
                if string:
                    fam_ids.append(rand_hog)
                    print(rand_hog)

        return fam_ids


def see_results(results):
    for k, v in results.items():
        print('{} queries found for {}'.format(len(v), k))


def get_hog_ids_from_results(results):

    query_event_hog_ids = {}

    for query, result in results.items():
        query_event_hog_ids[query] = [hashutils.result2fam(query)] + [hashutils.result2fam(r) for r in result]

    return query_event_hog_ids
"""
