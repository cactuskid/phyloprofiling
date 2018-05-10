import datasketch
import itertools

import hpputils


def hogid2fam(hog_id):
    fam = int(hog_id.split(':')[1])

    return fam

def fam2hogid(fam_id):
    """
    returns the hog fam id for any key of the lsh
    """
    hog_id = "HOG:" + (7-len(fam_id)) * '0' + str(fam_id)

    return hog_id

def result2hogid(result):
    """
    returns the hog fam id for any key of the lsh
    """
    hog_id = fam2hogid(result2fam(result))

    return hog_id

def result2fam(result):
    fam = str(result.split('-', 1)[0])

    return fam

def result2events(result):
    events = [event for event in result.split('-')[1:]]

    return events

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
    hashesDict = {}

    for event in eventdict:

        eventdict[event] = set(eventdict[event])
        minHash = datasketch.MinHash(num_perm=128)

        for element in eventdict[event]:
            minHash.update(element.encode())

        hashesDict[event] = minHash

    return hashesDict


def minhashes2leanminhashes(fam, minhashes, combination):

    lminhashesDict = {}

    for name, minhash in minhashes:
        lminHash = datasketch.LeanMinHash(minhash)
        lminHashName = str(fam) + '-' + name
        lminhashesDict[lminHashName] = lminHash

    if combination:
        for j in range(1, len(minhashes.keys())):
            for i in itertools.combinations(minhashes.keys(), j + 1):
                combName = str(fam)
                minHash = datasketch.MinHash(num_perm=128)
                for array in i:
                    combName += '-' + array
                    minHash.merge(minhashes[array])

                lminHash = datasketch.LeanMinHash(minHash)
                lminhashesDict[combName] = lminHash

    return lminhashesDict


def tree2hashes(fam, treemap, events, combination):

    eventdict = tree2eventdict(treemap)
    eventdict = {e: eventdict[e] for e in events}
    minhashes = eventdict2minhashes(eventdict)
    leanminhashes = minhashes2leanminhashes(fam, minhashes, combination)

    return leanminhashes


def fam2hashes(fam, dbObj, tree, replacement_dic, events, combination):
    """
    Turn the family into minhash object
    :param fam:
    :return:
    """
    treemap = hpputils.get_ham_treemap(fam, dbObj, tree, replacement_dic)
    hashes = tree2hashes(fam, treemap, events, combination)

    return hashes

