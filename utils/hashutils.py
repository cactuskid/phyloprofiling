import datasketch
import itertools
from scipy.sparse import csr_matrix


def hogid2fam(hog_id):
    """
    Get fam given hog id
    :param hog_id: hog id
    :return: fam
    """
    fam = int(hog_id.split(':')[1])

    return fam


def fam2hogid(fam_id):
    """
    Get hog id given fam
    :param fam_id: fam
    :return: hog id
    """
    hog_id = "HOG:" + (7-len(str(fam_id))) * '0' + str(fam_id)

    return hog_id


def result2hogid(result):
    """
    Get hog id given result
    :param result: result
    :return: hog id
    """
    hog_id = fam2hogid(result2fam(result))

    return hog_id


def result2fam(result):
    """
    Get fam given result
    :param result: result
    :return: fam
    """
    fam = int(result.split('-')[0])

    return fam


def result2events(result):
    """
    Get list of events given result
    :param result: result
    :return: list of events
    """
    events = [event for event in result.split('-')[1:]]

    return events


def tree2eventdict(treemap):
    """
    Get events dictionary from treemap object from pyham
    :param treemap: treemap object from pyham
    :return: dictionary of evolutionary events
    """
    eventdict = {'presence': [], 'gain': [], 'loss': [], 'duplication': []}
    for node in treemap.traverse():
        if not node.is_root():
            if node.nbr_genes > 0:
                eventdict['presence'].append('P' + node.name)
            if node.dupl > 0:
                eventdict['duplication'].append('D' + node.name)
            if node.lost > 0:
                eventdict['loss'].append('L' + node.name)
        else:
            eventdict['gain'].append('G' + node.name)
    return eventdict


def eventdict2minhashes(eventdict, nperm = 128):
    """
    Get minhashes from events dictionary
    :param eventdict: dictionary of evolutionary events
    :return: minHashes
    """
    hashes_dictionary = {}

    for event in eventdict:

        eventdict[event] = set(eventdict[event])
        minHash = datasketch.MinHash(num_perm=nperm)

        for element in eventdict[event]:
            minHash.update(element.encode())

        hashes_dictionary[event] = minHash

    return hashes_dictionary


def eventdict2WeightedMinhashes(eventdict, nperm = 128, treeweights =None ):
    """
    Get minhashes from events dictionary
    :param eventdict: dictionary of evolutionary events
    :return: minHashes
    """
	if wmg == None:
		wmg = WeightedMinHashGenerator()
	
	hashes_dictionary = {}

    for event in eventdict:
        eventdict[event] = set(eventdict[event])
        minHash = datasketch.MinHash(num_perm=nperm)

        for element in eventdict[event]:
            minHash.update(element.encode())

        hashes_dictionary[event] = minHash

    return hashes_dictionary



def minhashes2leanminhashes(fam, minhashes, combination=True ,  nperm = 128):
    """
    Get leanMinHashes from MinHashes
    :param fam: fam
    :param minhashes: list of minHashes
    :param combination: boolean, combination wanted or not
    :return: dictionary of leanMinHashes
    """
    lean_minhashes_dictionary = {}

    for name, minhash in minhashes.items():
        lean_minhash = datasketch.LeanMinHash(minhash)
        lean_minhash_name = str(fam) + '-' + name
        lean_minhashes_dictionary[lean_minhash_name] = lean_minhash

    if combination:
        for j in range(1, len(minhashes.keys())):
            for i in itertools.combinations(minhashes.keys(), j + 1):
                comb_name = str(fam)
                minHash = datasketch.MinHash(num_perm=nperm)
                for array in i:
                    comb_name += '-' + array
                    minHash.merge(minhashes[array])

                lean_minhash = datasketch.LeanMinHash(minHash)
                lean_minhashes_dictionary[comb_name] = lean_minhash

    return lean_minhashes_dictionary


def tree2hashes(fam, treemap, events, combination, nperm = 128):
    """
    Get hashes from tree
    :param fam: fam
    :param treemap: pyham treemap object
    :param events: list of events
    :param combination: boolean, combination wanted or not
    :return: dictionary containing two lists: one of minHashes and one of LeanMinHashes
    """
    if treemap is not None:
        event_dictionary = tree2eventdict(treemap)
        event_dictionary = {e: event_dictionary[e] for e in events}
        minhashes = eventdict2minhashes(event_dictionary,nperm)
        leanminhashes = minhashes2leanminhashes(fam, minhashes, combination, nperm)
        return {'hashes': minhashes, 'dict': leanminhashes}
    else:
        return None


def tree2hashes_from_row(row, events, combination,  nperm = 128):
    """
    Get hashes from tree
    :param row: tumple containing fam and treemap (pyham object)
    :param events: list of events
    :param combination: boolean, combination wanted or not
    :return: dictionary containing two lists: one of minHashes and one of LeanMinHashes
    """
    fam, treemap = row.tolist()

    return tree2hashes(fam, treemap, events, combination , nperm)


def tree2mat(treemap, taxa_index, verbose=False):
    """
    Turn each tree into a sparse matrix with 4 rows
    :param treemap: tree profile object from pyHam; contains biological events
    :param taxa_index: index of global taxonomic tree from OMA
    :param verbose: boolean, True for more info
    :return: profile_matrix : matrix of size numberOfBiologicalEvents times taxaIndex containing
    when the given biological event is present in the given species
    """
    if treemap is not None:
        hog_matrix = csr_matrix((1, 4 * len(taxa_index)))
        column_dict = {'presence': 0, 'gain': 1, 'loss': 2, 'duplication': 3}

        for node in treemap.traverse():
            # traverse() returns an iterator to traverse the tree structure
            # strategy:"levelorder" by default; nodes are visited in order from root to leaves
            # it return treeNode instances
            if not node.is_root():
                # for presence, loss, and duplication, set 1 if > 0
                if node.nbr_genes > 0:
                    pad = column_dict['presence'] * len(taxa_index)
                    hog_matrix[0, pad + taxa_index[node.name]] = node.nbr_genes
                if node.lost > 0:
                    pad = column_dict['loss'] * len(taxa_index)

                    hog_matrix[0, pad + taxa_index[node.name]] = 1
                if node.dupl > 0:
                    pad = column_dict['duplication'] * len(taxa_index)
                    hog_matrix[0, pad + taxa_index[node.name]] = 1
            else:
                # gain is only for root; impossible to "gain" a gene several times
                pad = column_dict['gain'] * len(taxa_index)
                hog_matrix[0, pad + taxa_index[node.name]] = 1

        if verbose:
            print(hog_matrix)

        return hog_matrix

    else:
        return csr_matrix((1, 4 * len(taxa_index)))


def fam2hash_hdf5(fam, hdf5, events=['duplication', 'gain', 'loss', 'presence']):
    """
    Get the minhash corresponding to the given hog id number
    :param fam: hog id number
    :param hdf5: hdf5 file containing hashes
    :param events: list of events the hashes are build on; default: all four events
    :return: list of hashes for the given fam and events
    """
    minhash1 = None
    for event in events:
        query_minhash = datasketch.MinHash(num_perm=128)

        try:
            query_minhash.hashvalues = hdf5[event][fam, :]
        except:
            print('bug fam2hash_hdf5')
            print(fam)
            print(event)
            print(hdf5[event].shape)

        query_minhash.seed = 1

        if minhash1 is None:
            minhash1 = datasketch.MinHash(seed=query_minhash.seed, hashvalues=query_minhash.hashvalues)
        else:
            minhash2 = datasketch.MinHash(seed=query_minhash.seed, hashvalues=query_minhash.hashvalues)
            minhash1.merge(minhash2)

    return minhash1
