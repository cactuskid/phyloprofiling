
# use dask and pyham to create a big sparse matrix for each type of evolutionary event
import ete3
from scipy.sparse import csr_matrix,find , vstack
from datasketch import MinHash, LeanMinHash
import numpy as np
import config
import datasketch

def generateTaxaIndex(species_tree):
    '''
    Generates an index for the global taxonomic tree for all OMA
    Args:
        species_tree: species tree in newick format
    Returns:
        taxaIndex: dictionary key: node name (species name); value: index
        taxaIndexReverse: dictonary key: index: value: species name
    '''
    t = ete3.Tree(species_tree, format=1)
    taxaIndex = {}
    taxaIndexReverse = {}
    for i,node in enumerate(t.traverse()):
        taxaIndexReverse[i] = node.name
        taxaIndex[node.name] = i
    return taxaIndex, taxaIndexReverse


def FamList2RowsOnTheFly(listfam, dbObj, tree, dic):
    taxaIndex, taxaIndexReverse = generateTaxaIndex(tree)
    rows = []

	for fam in listfam:
		# get treemap
		treemap_fam = pyhamPipeline.get_hamTree(fam, dbObj, tree, dic)
		# get mat
		rows.append(Tree2mat(treemap_fam, taxaIndex))
	stackRows = vstack(rows)
	return stackRows


def FamList2RowsH5(h5file, listfam):
    pass


def jaccard_cutoff(fams, scores, cutoff):
    return fams[ np.where(scores > cutoff)]

# TODO debug this
def get_hash_hog_id(fam, h5hashes, events=['duplication', 'gain', 'loss', 'presence']):
    # get hash of desired events
    print('test')
    minhash1 = None
    for event in events:
        query_minhash = LeanMinHash(MinHash(num_perm=128))
        query_minhash.hashvalues = h5hashes[event][fam, :]
        query_minhash.seed = 1
        # np.get_buffer doesnt work
        #buf = bytearray(h5hashes[event][fam,:])
        #buf = np.get_buffer(h5hashes[event][fam,:])
        #buf = np.frombuffer(h5hashes[event][fam,:])
        if minhash1 is None:
            #query_minhash = LeanMinHash.deserialize(buf)
            minhash1 = MinHash(seed=query_minhash.seed, hashvalues=query_minhash.hashvalues)
        else:
            #query_minhash = LeanMinHash.deserialize(buf)
            minhash2 = MinHash(seed=query_minhash.seed, hashvalues=query_minhash.hashvalues)
            minhash1.merge(minhash2)
    return minhash1

def test_gethoghash():
    with File( config.datadir + 'hashes.h5', 'r') as h5hashes:
        entries = range(100,200)
        for i in entries:
            print(get_hash_hog_id(i, h5hashes))


def get_hashdict(fams, h5hashes, events = ['duplication', 'gain', 'loss', 'presence']):
    hashdict = {}
    for fam in fams:
        hashdict[fam] = get_hash_hog_id( fam, h5hashes, events)
    return hashdict

def jaccard_rank(query,hashdict):
	hashlist = np.asarray(list(hashdict.values()))
	fams = np.asarray(list(hashdict.keys()))
	scores = [ queryhash.jaccard(resulthash) for resulthash in hashlist ]
	index = np.argsort(scores)
	for array in [hashlist, fams , index]:
		array = array[index]
	return fams,scores


def Tree2Hashes(treemap, fam=None, LSH=None , l = None):

    #turn each tree into a minhash object
    #serialize and store as array
    eventdict = { 'presence':[] , 'gain':[] , 'loss':[] , 'duplication':[]}
    for node in treemap.traverse():
    # traverse() returns an iterator to traverse the tree structure
    # strategy:"levelorder" by default; nodes are visited in order from root to leaves
    # it return treeNode instances
        if not node.is_root():
            if node.nbr_genes >0:
                eventdict['presence'].append('P'+node.name)
            if node.dupl > 0:
                eventdict['duplication'].append('D'+node.name)
            if node.lost > 0:
                eventdict['loss'].append('L'+node.name)
        else:
            eventdict['gain'].append('G'+node.name)

    hashes = {}

    hashesDict = {}

    lminHashDict = {}

    for array in eventdict:
        eventdict[array] = set(eventdict[array])

        minHash = datasketch.MinHash(num_perm=128)

        for element in eventdict[array]:

            minHash.update(element.encode())


        hashesDict[array] = minHash

        lminHash = datasketch.LeanMinHash(minHash)

        lminHashName = str(fam)+'-'+array

        lminHashDict[lminHashName] = lminHash
        if LSH is not None:
            LSH.insert(lminHashName, lminHash)


        buf = bytearray(lminHash.bytesize())
        lminHash.serialize(buf)

        hashes[array] = buf


    #hashmat = np.hstack(buffers)

    for j in range(1,len(eventdict.keys())):
        for i in itertools.combinations(eventdict.keys(), j+1):
            combName = str(fam)
            minHash = datasketch.MinHash(num_perm=128)
            for array in i:
                combName += '-'+array
                minHash.merge(hashesDict[array])

            lminHashDict[combName] = minHash
            lminHash = datasketch.LeanMinHash(minHash)
        if LSH:
            LSH.insert( combName , lminHash)

    return hashes , lminHashDict

def DFTree2Hashes(row):
    #turn each tree into a minhash object
    #serialize and store as array

    fam, treemap = row.tolist()

    eventdict = { 'presence':[] , 'gain':[] , 'loss':[] , 'duplication':[]}
    if treemap is not None:
        for node in treemap.traverse():
        # traverse() returns an iterator to traverse the tree structure
        # strategy:"levelorder" by default; nodes are visited in order from root to leaves
        # it return treeNode instances
            if not node.is_root():
                if node.nbr_genes >0:
                    eventdict['presence'].append('P'+node.name)
                if node.dupl > 0:
                    eventdict['duplication'].append('D'+node.name)
                if node.lost > 0:
                    eventdict['loss'].append('L'+node.name)
            else:
                eventdict['gain'].append('G'+node.name)

        hashes = {}

        hashesDict = {}

        lminHashDict = {}

        for array in eventdict:
            eventdict[array] = set(eventdict[array])

            minHash = datasketch.MinHash(num_perm=128)

            for element in eventdict[array]:

                minHash.update(element.encode())


            hashesDict[array] = minHash

            lminHash = datasketch.LeanMinHash(minHash)

            lminHashName = str(fam)+'-'+array

            lminHashDict[lminHashName] = lminHash

            #lminHash.bytesize() --> computes the bitesize
            #buf = bytearray(lminHash.bytesize())
            #lminHash.serialize(buf)

            #hashes[array] = buf

            hashes[array] = lminHash.hashvalues


			buf = bytearray(lminHash.bytesize())
			lminHash.serialize(buf)

			hashes[array] = buf
        #hashmat = np.hstack(buffers)

        for j in range(1,len(eventdict.keys())):
            for i in itertools.combinations(eventdict.keys(), j+1):
                combName = str(fam)
                minHash = datasketch.MinHash(num_perm=128)
                for array in i:
                    combName += '-'+array
                    minHash.merge(hashesDict[array])

                lminHashDict[combName] = minHash
                lminHash = datasketch.LeanMinHash(minHash)

        return {'hashes':hashes , 'dict': lminHashDict}
    else:
        return None

def Tree2mat( treemap, taxaIndex, verbose = False):
    '''
    Turn each tree into a sparse matrix with 4 rows.

    Args:
        treemap : tree profile object from pyHam; contains biological events
        taxaIndex : index of gloval taxonomic tree from OMA

    Returns :
        profile_matrix : matrix of size numberOfBiologicalEvents times taxaIndex containing when the given biological event is present in the given species
    '''
    #use partials to configure the taxa index
    if treemap is not None:
        hogmat = csr_matrix( (1 , 4*len(taxaIndex)  ) )
        columndict={ 'presence':0
         , 'gain': 1
         , 'loss': 2
         , 'duplication':3
        }


        for node in treemap.traverse():
        # traverse() returns an iterator to traverse the tree structure
        # strategy:"levelorder" by default; nodes are visited in order from root to leaves
        # it return treeNode instances
            if not node.is_root():
                # for presence, loss, and duplication, set 1 if > 0
                if node.nbr_genes > 0:
                    pad = columndict['presence'] *len(taxaIndex)
                    hogmat[ 0 , pad + taxaIndex[node.name]] = node.nbr_genes
                if node.lost > 0:
                    pad = columndict['loss']*len(taxaIndex)

                    hogmat[0 , pad + taxaIndex[node.name]] = 1
                if node.dupl > 0:
                    pad = columndict['duplication']*len(taxaIndex)
                    hogmat[0 , pad + taxaIndex[node.name]] = 1
            else:
                # gain is only for root; impossible to "gain" a gene several times
                pad = columndict['gain']*len(taxaIndex)
                hogmat[0 , pad + taxaIndex[node.name]] = 1

        if verbose == True:
            print(hogmat)

        return hogmat

    else:
        return csr_matrix( (1 , 4*len(taxaIndex)  ) )
