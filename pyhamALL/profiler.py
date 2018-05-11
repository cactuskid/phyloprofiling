import pickle
import tables
import pandas as pd
import numpy as np
import multiprocessing as mp
import itertools

from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import TermCounts
from pyoma.browser import db
import pyham

import validation.phyloValidationGoTerm as validationGoTerm
import hpputils
import hashutils



class Profiler:

    def __init__(self, h5_oma_path):
        # load all the files (oma database, go dag file, association gaf file=
        self.h5OMA = tables.open_file(h5_oma_path, mode='r')
        self.dbObj = db.Database(self.h5OMA)
        self.omaIdObj = db.OmaIdMapper(self.dbObj)
        self.replacement_dic, self.tree = hpputils.create_species_tree(self.h5OMA, self.omaIdObj)

    def go_benchmarking_init(self, obo_file_path, gaf_file_path):
        self.go = obo_parser.GODag(obo_file_path)
        self.associations = read_gaf(gaf_file_path)
        # Get the counts of each GO term.
        self.termcounts = TermCounts(self.go, self.associations)
        self.goTermAnalysis = validationGoTerm.SemanticSimilarityAnalysis(self.go, self.h5OMA, self.termcounts)

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
        query_hashes = hashutils.fam2hashes(fam_id, self.dbOjb, self.tree, self.replacement_dic, events, combination)

        resultDict = {}
        for name, hashvalue in query_hashes:
            # lsh.query returns a list of hash (fam-event1-event2- etc)
            resultDict[name] = self.lsh.query(hashvalue)

        return resultDict

    def results(self, hog_id=None, fam_id=None, events=['duplication', 'gain', 'loss', 'presence'], combination=True, scores=False):

        results = self.hog_query(hog_id=hog_id, fam_id=fam_id, events=events, combination=combination)

        if not scores:
            results_df = pd.DataFrame.from_dict(results)
            return results_df

        result_dict = {}

        # semantic scores
        for query, result in results:
            related_hog_list = list(set([query] + [r for r in result]))
            events_combo = hashutils.result2events(query)
            result_dict = self.compute_semantic_distance(related_hog_list, result_dict, events_combo)

        # get jaccard score, here its possible to use up
        result_dict = self.compute_jaccard_distance(related_hog_list, result_dict, events_combo)

        return result_dict

    def compute_semantic_distance(self, hogs_list, result_dict, events_combo):

        for i, hog_event_1 in enumerate(hogs_list):
            for j, hog_event_2 in enumerate(hogs_list):
                hog1 = hashutils.result2hogid(hog_event_1)
                hog2 = hashutils.result2hogid(hog_event_2)

                # get semantic dist between two hog (ids) and save it in matrix and in big dict
                semantic_dist = self.goTermAnalysis.semantic_similarity_score(hog1, hog2)
                result_dict[(hog1, hog2, events_combo)]['semantic_distance'] = semantic_dist
                # also save the go terms in big dict
                result_dict[(hog1, hog2, events_combo)]['go term hog 1'] = self.goTermAnalysis.get_go_terms(hog1)
                result_dict[(hog1, hog2, events_combo)]['go term hog 2'] = self.goTermAnalysis.get_go_terms(hog2)


        return result_dict

    def jaccard_mp(self, hash1, hash2):
        """
        computes jaccard score for two hashes
        returns score and index
        """
        i, h1 = hash1
        j, h2 = hash2

        return i, j, h1.jaccard(h2)

    def compute_jaccard_distance(self, hogs_list, result_dict, events_combo):
        # prepare data for mp stuff
        pool = mp.Pool()

        hashes = {result: hashutils.result2hash(result) for result in hogs_list}

        hashcompare = [(h1, h2) for h1, h2 in itertools.combinations(enumerate(hashes), 2)]
        # compute jaccard
        jaccard_results = pool.map_async(hashcompare, jaccard_mp).get()
        # feed matrix and big dict with results
        for jac_res in jaccard_results:
            i, j, jac_dist = jac_res
            hog1 = hashutils.result2hogid(hogs_list[i])
            hog2 = hashutils.result2hogid(hogs_list[j])
            result_dict[(hog1, hog2, events_combo)]['jaccard_distance'] = jac_dist

        return result_dict

    def results_save(self, hog_id=None, fam_id=None, events=['duplication', 'gain', 'loss', 'presence'],
                    combination=True, scores=False, path_to_save=None):
        df_results = self.results(hog_id=hog_id, fam_id=fam_id, events=events, combination=combination, scores=scores)

        if path_to_save is not None:
            df_results.to_csv(path_to_save, sep='\t')