
import _pickle as pickle
import pandas as pd
import h5py
import itertools
import ujson as json
import random
from scipy.sparse import csr_matrix
from tables import *
from pyoma.browser import db
import numpy as np
import random
np.random.seed(0)
random.seed(0)
import ete3
from datasketch import WeightedMinHashGenerator
from validation import validation_semantic_similarity
from utils import hashutils,  config_utils , pyhamutils , files_utils
from time import time

class Profiler:


    def __init__(self,lshforestpath = None, hashes_h5=None, mat_path= None, omapath = None):
        #use the lsh forest or the lsh

        """
        A profiler object allows the user to query the LSH with HOGs and get a list of result HOGs back

        """
        print('loading lsh')
        with open(lshforestpath, 'rb') as lshpickle:
            self.lshobj = pickle.loads(lshpickle.read())
            print('indexing lsh')
            self.lshobj.index()
        self.hashes_h5 = h5py.File(hashes_h5, mode='r')
        print('DONE')

        if mat_path:

            ## TODO: change this to read hdf5
            #profile_matrix_file = open(profile_matrix_path, 'rb')
            #profile_matrix_unpickled = pickle.Unpickler(profile_matrix_file)
            #self.profile_matrix = profile_matrix_unpickled.load()
            pass

        if omapath:
            #open oma db object
            ## TODO: unfinished
            #open up taxa Index
            self.taxtree = ete3.Tree.phylotree('./mastertree.nwk')
            self.taxaIndex = { n.name:i for i,n in enumerate(self.taxtree.traverse()) }

            h5_oma = open_file(config_utils.omadir + 'OmaServer.h5', mode="r")
            self.db_obj = db.Database(h5_oma)

            #open up master tree
            self.HAM_PIPELINE = functools.partial(pyhamutils.get_ham_treemap_from_row, tree=self.tree_string )
            self.HASH_PIPELINE = functools.partial(hashutils.row2hash , taxaIndex=self.taxaIndex  , treeweights=self.treeweights , wmg=wmg )
            self.READ_ORTHO = functools.partial(pyhamutils.get_orthoxml, db_obj=self.db_obj)


            #setup function to generate a profile on the fly

    def return_profile_OTF(self, fam):
        ## TODO: unfinished
        ortho_fam = self.READ_ORTHO(fam)
        pyham_tree = self.HAM_PIPELINE([fam, ortho_fam])
        losses = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.lost and n.name in self.taxaIndex  ]
        dupl = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.dupl  and n.name in self.taxaIndex  ]
        presence = [ self.taxaIndex[n.name]  for n in tp.traverse() if n.nbr_genes > 0  and n.name in self.taxaIndex  ]
        return mat

    def return_profile_mat_OTF(self , hogs):
        return np.vstack([ self.return_profile_OTF(hog) for hog in hogs])


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

        results = self.lshobj.query(query_hash, k)
        return results

    def hog_query_OMA(self,hog_id=None, fam_id=None , k = 100 ):
        #construct a profile on the fly
        #rand seed values need to be identical between the construction of the lsh DB and the use of this function
        """
        Untested, Given a hog_id or a fam_id as a query, returns a dictionary containing the results of the LSH.
        Generates the tree profile and hashes on the fly
        :param hog_id: query hog id
        :param fam_id: query fam id
        :return: list containing the results of the LSH for the given query
        """
        if hog_id is not None:
            fam_id = hashutils.hogid2fam(hog_id)
        ortho = self.lshobj.READ_ORTHO(fam)
        tp = self.lshobj.HAM_PIPELINE((fam, ortho))
        hash = self.lshobj.HASH_PIPELINE((fam, tp))
        results = self.lshobj.query(query_hash, k)
        return results

    def pull_hashes(self , hoglist):
        """
        Given a list of hog_ids , returns a dictionary containing their hashes.
        This uses the hdf5 file to get the hashvalues
        :param hog_id: query hog id
        :param fam_id: query fam id
        :return: a dict containing the hash values of the hogs in hoglist
        """
        return { hog:hashutils.fam2hash_hdf5(hog, self.hashes_h5 ) for hog in hoglist}

    def pull_go(self,hoglist):
        pass

    def pull_matrows(self,fams):
        """
        given a list of fams return the submatrix containing their profiles

        :return:fams sorted, sparse mat
        """
        return self.profile_matrix[np.asarray(fams),:]
    @staticmethod
    def sort_hashes(query_hash,hashes):
        """
        Given a dict of hogs:hashes, returns a sorted array of hogs and jaccard distances relative to query hog.
        :param query hash: weighted minhash of the query
        :param hashes: a dict of hogs:hashes
        :return: sortedhogs, jaccard
        """
        #sort the hashes by their jaccard relative to query hash
        jaccard=[ query_hash.jaccard(hashes[hog]) for hog in hashes]
        index = np.argsort(jaccard)
        sortedhogs = np.asarry(list(hashes.keys()))[index]
        jaccard= jaccard[index]
        return sortedhogs, jaccard


    @staticmethod
    def allvall_hashes(hashes):
        """
        Given a dict of hogs:hashes, returns generate an all v all jaccard distance matrix.
        :param hashes: a dict of hogs:hashes
        :return: hashmat
        """
        #generate an all v all jaccard distance matrix
        hashmat = np.zeros((len(hashes),len(hashes)))
        for i , hog1 in enumerate(hashes):
            for j, hog2 in enumerate(hahes):
                hashmat[i,j]= hashes[hog1].jaccard(hashes[hog2])
        return hashmat

    @staticmethod
    def hog_v_hog(hog1,hog2):
        """
        give two hogs returns jaccard distance.
        :param hog1 , hog2: str hog id
        :return: jaccard score
        """
        #generate an all v all jaccard distance matrix
        hashes = self.pull_hashes([hog1,hog2])
        hashes = list(hashes.values())

        return hashes[0].jaccard(hashes[1])

    def allvall_nx(G,hashes,thresh =None):
        """
        Given a dict of hogs:hashes, returns generate an all v all jaccard distance matrix.
        :param hashes: a dict of hogs:hashes
        :return: hashmat
        """
        #generate an all v all jaccard distance matrix
        hashmat = np.zeros((len(hashes),len(hashes)))

        for i , hog1 in enumerate(hashes):
            for j, hog2 in enumerate(hashes):
                hashmat[i,j]= hashes[hog1].jaccard(hashes[hog2])
        return hashmat

    def iternetwork(seedHOG):
        pass

    def rank_hashes(query_hash,hashes):
        jaccard = []
        sorted = []
        scores = {}
        hogsRanked = np.asarray(list(hashes.keys()))
        for i, hog in enumerate(hashes):
            score = query_hash.jaccard(hashes[hog])
            jaccard.append( score)
            scores[hog] = score
        hogsRanked = list( hogsRanked[ np.argsort(jaccard) ] )
        jaccard = np.sort(jaccard)
        return hogsRanked, jaccard


    def get_vpairs(fam):

        """
        get pairwise distance matrix of OMA all v all
        #not finished
        :param fam: an oma fam
        :return sparesemat: a mat with all taxa in Oma with nonzero entries where this protein is found
        :return densemat: a mat with the taxa covered by the fam
        """
        taxa = self.db_obj.hog_levels_of_fam(fam)
        subtaxindex = { taxon:i for i,taxon in enumerate(taxa)}
        prots = self.db_obj.hog_members_from_hog_id(fam,  'LUCA')
        for prot in prots:
            taxon = prot.ncbi_taxon_id()
            pairs = self.db_obj.get_vpairs(prot)
            for EntryNr1, EntryNr2, RelType , score , distance in list(pairs):
                pass
        return sparsemat , densemat


    def get_submatrix_form_results(self, results):
        res_mat_list = []
        for query, result in results.items():
            res_mat = csr_matrix((len(result), self.profile_matrix.shape[1]))
            for i, r in enumerate(result):
                res_mat[i, :] = self.profile_matrix[r, :]
            res_mat_list.append(res_mat)
        final = np.vstack(res_mat_list)
        return res_mat_list
