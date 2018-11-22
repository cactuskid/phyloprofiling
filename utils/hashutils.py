

import datasketch
import itertools
from scipy.sparse import lil_matrix
import ete3
import copy
import math
import numpy as np
import pandas as pd


def generate_treeweights( mastertree, taxaIndex , taxfilter, taxmask , lambdadict, start, exp = True):
    #weighing function for tax level, masking levels etc
    """
    Generate the weights of each taxonomic level to be applied during the
    constructin of weighted minhashes

    :param mastertree: full corrected ncbi taxonomy
    :param taxaIndex: dict mapping taxa to columns
    :param taxfilter: list of branches to delete
    :param taxmask: if this is not NONE taxmask, the DB is constructed with this subtree
    :param lambdadict: parameters for weight functions
    :param start: parameters for weight functions
    :return: weights: a vector of weights for each tax level
    """

    weights = { type: np.zeros((len(taxaIndex),1)) for type in ['presence', 'loss', 'dup']}
    for node in mastertree.traverse():
        node.add_feature('degree', 1 )
    for node in mastertree.iter_descendants():
        for d in node.iter_descendants():
            d.degree+=1
    for event in weights:
        newtree = copy.deepcopy(mastertree)
        for n in newtree.traverse():
            if taxmask:
                if str(n.name) == str(taxmask):
                    newtree = n
                    break
            if taxfilter:
                if n.name in taxfilter:
                    #set weight for descendants of n to 0
                    n.delete()
        for n in newtree.traverse():
            #exponential decay of initial weigh over node degree
            if exp == True:
                #exponential
                weights[event][taxaIndex[n.name]] = start[event]*math.exp(n.degree*lambdadict[event])
            else :
                #linear
                weights[event][taxaIndex[n.name]] = start[event] + n.degree*lambdadict[event]
    return weights

def hash_tree(tp , taxaIndex , treeweights , wmg):
    """
    Generate the weights of each taxonomic level to be applied during the
    constructin of weighted minhashes

    :param tp: a pyham tree profile
    :param taxaIndex: dict mapping taxa to columns
    :param treeweights: a vector of weights for each tax levels
    :param wmg: Datasketch weighted minhash generator
    :return hog_matrix: a vector of weights for each tax level
    :return weighted_hash: a weighted minhash of a HOG

    """
    #convert a tree profile to a weighted minhash
    losses = [ taxaIndex[n.name]  for n in tp.traverse() if n.lost and n.name in taxaIndex  ]

    dupl = [ taxaIndex[n.name]  for n in tp.traverse() if n.dupl  and n.name in taxaIndex  ]

    presence = [ taxaIndex[n.name]  for n in tp.traverse() if n.nbr_genes > 0  and n.name in taxaIndex  ]

    indices = dict(zip (['presence', 'loss', 'dup'],[presence,losses,dupl] ) )

    hog_matrix = lil_matrix((1, 3*len(taxaIndex)))
    hogsum = 0
    for i,event in enumerate(indices):
        if len(indices[event])>0:
            taxindex = np.asarray(indices[event])
            hogindex = np.asarray(indices[event])+i*len(taxaIndex)
            #print(index)
            #print(treeweights[event][index].shape)
            #print(hog_matrix[:,index].shape)
            #assign the
            hog_matrix[:,hogindex] = treeweights[event][taxindex].ravel()
            hogsum+=np.sum(treeweights[event][taxindex])

    #normalize total...
    if hogsum > 0:
        weighted_hash = wmg.minhash(list(hog_matrix.todense().flat))
        return  hog_matrix,weighted_hash
    else:
        return None, None

def row2hash(row , taxaIndex , treeweights , wmg):
    """
    turn a dataframe row with an orthoxml file to hash and matrix row
    :param row: lsh builder dataframe row
    :param taxaIndex: dict mapping taxa to columnsfam
    :param treeweights: a vector of weights for each tax levels
    :param wmg: Datasketch weighted minhash generator
    :return: hog_matrix: a vector of weights for each tax level
    :return: weighted_hash: a weighted minhash of a HOG
    """
    #convert a dataframe row to a weighted minhash
    fam, treemap = row.tolist()
    hog_matrix,weighted_hash = hash_tree(treemap , taxaIndex , treeweights , wmg)
    return [weighted_hash,hog_matrix]


def fam2hash_hdf5(fam,  hdf5, dataset = None, nsamples = 128  ):
    #read the stored hash values and return a weighted minhash
    """
    Read the stored hash values and return a weighted minhash
    :param fam: hog id
    :param hdf5: h5py object of the hashvalues
    :param dataset: which dataset to use when constructing the hash
    :return: minhash1: the weighted hash of your HOG
    """
    if dataset is None:
        #use first dataset by default
        dataset = list(hdf5.keys())[0]
    hashvalues = np.asarray(hdf5[dataset][fam, :].reshape(nsamples,2 ))
    hashvalues = hashvalues.astype('int64')
    minhash1 = datasketch.WeightedMinHash( seed = 1, hashvalues=hashvalues)

    return minhash1


def hogid2fam(hog_id):
    """
    Get fam given hog id
    :param hog_id: hog id
    :return: fam
    """
    if ':' in hog_id:
        fam = int(hog_id.split(':')[1])
    else:
        fam = int(hog_id)
    return fam


def fam2hogid(fam_id):
    """
    Get hog id given fam
    :param fam_id: fam
    :return: hog id
    """
    hog_id = "HOG:" + (7-len(str(fam_id))) * '0' + str(fam_id)

    return hog_id
