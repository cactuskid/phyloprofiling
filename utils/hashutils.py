import datasketch
import itertools
from scipy.sparse import csr_matrix

from utils import files_utils


def hogid2fam(hog_id):
    fam = int(hog_id.split(':')[1])

    return fam


def fam2hogid(fam_id):
    """
    returns the hog fam id for any key of the lsh
    """

    hog_id = "HOG:" + (7-len(str(fam_id))) * '0' + str(fam_id)

    return hog_id


def result2hogid(result):
    """
    returns the hog fam id for any key of the lsh
    """
    hog_id = fam2hogid(result2fam(result))

    return hog_id


def result2fam(result):
    fam = int(result.split('-', 1)[1])

    return fam


def result2events(result):
    events = [event for event in result.split('-')[:-1]]

    return events


def result2hash(result, tree, replacement_dic, db_obj):
    events = result2events(result)
    fam = result2fam(result)

    treemap = files_utils.get_ham_treemap(fam, db_obj, tree, replacement_dic)
    eventdict = tree2eventdict(treemap)
    eventdict = {e: eventdict[e] for e in events}

    minhashes = eventdict2minhashes(eventdict)

    hash_value = combine_minhashes(minhashes)

    return hash_value


def combine_minhashes(hashes):
    minhash = datasketch.MinHash(num_perm=128)
    for name, hashvalue in hashes:
        minhash.merge(hashvalue)

    return minhash


def tree2eventdict(treemap):
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


def eventdict2minhashes(eventdict):
    hashes_dictionary = {}

    for event in eventdict:

        eventdict[event] = set(eventdict[event])
        minHash = datasketch.MinHash(num_perm=128)

        for element in eventdict[event]:
            minHash.update(element.encode())

        hashes_dictionary[event] = minHash

    return hashes_dictionary


def minhashes2leanminhashes(fam, minhashes, combination):

    lean_minhashes_dictionary = {}

    for name, minhash in minhashes:
        lean_minhash = datasketch.LeanMinHash(minhash)
        lean_minhash_name = str(fam) + '-' + name
        lean_minhashes_dictionary[lean_minhash_name] = lean_minhash

    if combination:
        for j in range(1, len(minhashes.keys())):
            for i in itertools.combinations(minhashes.keys(), j + 1):
                combName = str(fam)
                minHash = datasketch.MinHash(num_perm=128)
                for array in i:
                    combName += '-' + array
                    minHash.merge(minhashes[array])

                lean_minhash = datasketch.LeanMinHash(minHash)
                lean_minhashes_dictionary[combName] = lean_minhash

    return lean_minhashes_dictionary


def tree2hashes(fam, treemap, events, combination):

    if treemap is not None:
        event_dictionary = tree2eventdict(treemap)
        event_dictionary = {e: event_dictionary[e] for e in events}
        minhashes = eventdict2minhashes(event_dictionary)
        leanminhashes = minhashes2leanminhashes(fam, minhashes, combination)

        return leanminhashes

    else:
        return None


def tree2hashes_from_row(row, events, combination):
    fam, treemap = row.tolist()

    return tree2hashes(fam, treemap, events, combination)


def fam2hashes(fam, db_obj, tree, replacement_dic, events, combination):
    """
    Turn the family into minhash object
    :param fam: hog family id
    :param db_obj: database object
    :param tree: tree
    :param replacement_dic: replacement dictionary
    :param events: list of evolutionary events
    :param combination: boolean: True for combination of events
    :return: hashes corresponding to this tree
    """
    treemap = files_utils.get_ham_treemap(fam, db_obj, tree, replacement_dic)
    hashes = tree2hashes(fam, treemap, events, combination)

    return hashes


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