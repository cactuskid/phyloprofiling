
import _pickle as pickle
import pandas as pd
import h5py
import itertools
import ujson as json
import random
from scipy.sparse import csr_matrix

import lshbuilder

from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import TermCounts
from pyoma.browser import db

import numpy as np
import random

np.random.seed(0)
random.seed(0)


from datasketch import WeightedMinHashGenerator
from validation import validation_semantic_similarity
from utils import hashutils, string_utils , config_utils

from time import time
class Profiler:

    def __init__(self,lshforestpath = None, hashes_h5=None, mat_path= None,  lsh_builder_path = None, unimap_path = None ,string_data_path = None , GO= None):
        #use the lsh forest or the lsh
        if lsh_builder_path:
            with open( lsh_builder_path, mode='rb', buffering=None) as lshin:
                lsh_builder = pickle.loads(lshin.read())
            self.lsh_builder = lsh_builder

            self.lshpath = self.saving_path + 'newlsh.pkl'
            self.lshforestpath = self.saving_path + 'newlshforest.pkl'
            print('loading lsh')
            with open(lsh_builder.lshforestpath, 'rb') as lshpickle:
                self.lshobj = pickle.loads(lshpickle.read())
                print('indexing lsh')
                self.lshobj.index()

            self.hashes_h5 = h5py.File(lsh_builder.hashes_h5, mode='r')
            print('DONE')
        else:

            with open(lshforestpath, 'rb') as lshpickle:
                self.lshobj = pickle.loads(lshpickle.read())
                print('indexing lsh')
                self.lshobj.index()
            self.hashes_h5 = h5py.File(hashes_h5, mode='r')
            print('DONE')

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
        if GO :
            self.go = obo_parser.GODag(obo_file_path)
            self.associations = read_gaf(gaf_file_path)

            self.term_counts = TermCounts(self.go, self.associations)
            self.goTermAnalysis = validation_semantic_similarity.Validation_semantic_similarity(self.go,
                                                                                                self.term_counts,
                                                                                                self.go_terms_hdf5)


    def hog_query(self, hog_id=None, fam_id=None , k = 100 ):
        """
        Given a hog_id or a fam_id as a query, returns a dictionary containing the results of the LSH.
        :param hog_id: query hog id
        :param fam_id: query fam id
        :return: list containing the results of the LSH for the given query
        """
        if hog_id is not None:
            fam_id = hashutils.hogid2fam(hog_id)
        query_hash = hashutils.fam2hash_hdf5(fam_id, self.hashes_h5 , nsamples=  128 )
        #print(query_hash)

        results = self.lshobj.query(query_hash, k)
        return results



    def hog_query_OMA(hog_id=None, fam_id=None , k = 100 ):
        #construct a profile on the fly
        if hog_id is not None:
            fam_id = hashutils.hogid2fam(hog_id)

        ortho = self.lshobj.READ_ORTHO(fam)
        tp = self.lshobj.HAM_PIPELINE((fam, ortho))
        hash = self.lshobj.HASH_PIPELINE((fam, tp))

        results = self.lshobj.query(query_hash, k)

        return results

    def hog_query_manual(famtree , k = 100 ):
        #construct a profile on the fly
        #feed the hash pipeline a tree with events.
        #hacky way of getting
        hash = self.lshobj.HASH_PIPELINE(('1', tp))
        results = self.lshobj.query(query_hash, k)
        return results

    def pull_hashes(self , hoglist):
        return { hog:hashutils.fam2hash_hdf5(hog, self.hashes_h5 ) for hog in hoglist}

    def pull_mapping(self, hoglist):
        #grab the crossrefs for a list of hogs
        return { hog:{ dataset: json.loads(unimap_h5[dataset][hog]) for dataset in unimap_h5 }  for hog in hoglist }

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
        taxa = self.db_obj.hog_levels_of_fam(fam)
        subtaxindex = { taxon:i for i,taxon in enumerate(taxa)}
        prots = self.db_obj.hog_members_from_hog_id(fam,  'LUCA')
        for prot in prots:
            taxon = prot.ncbi_taxon_id()
            pairs = self.db_obj.get_vpairs(prot)
            for EntryNr1, EntryNr2, RelType , score , distance in list(pairs):
                pass
        return sparsemat , densemat


    def compute_semantic_distance(self, hog_1, hog_2):
        semantic_dist = self.goTermAnalysis.semantic_similarity_score(hog_1, hog_2)
        return semantic_dist

    def get_hogs_with_annotations(self):
        print('getting hogs')
        hogs_with_annotations = []
        for fam, goterms in enumerate(self.hogs2goterms):
            if goterms:
                if json.loads(goterms):
                    hogs_with_annotations.append(fam)
        print(len(hogs_with_annotations))
        return hogs_with_annotations


    def get_submatrix_form_results(self, results):

        res_mat_list = []

        for query, result in results.items():
            res_mat = csr_matrix((len(result), self.profile_matrix.shape[1]))
            for i, r in enumerate(result):
                res_mat[i, :] = self.profile_matrix[r, :]

            res_mat_list.append(res_mat)

        return res_mat_list
