import itertools

import redis

from utils import hashutils


def connect2IDmap():
    r1 = redis.StrictRedis(host='10.0.63.33', port=6379, db=0)
    return r1


def connect2Stringmap():
    r2 = redis.StrictRedis(host='10.0.63.33', port=6379, db=1)
    return r2


def fam2stringID(dbobj, fam, r):
    members = dbobj.iter_members_of_hog_id(hashutils.fam2hogid(fam))
    allstring = [r.get(m.canonicalid) for m in members]
    return allstring


def HOGvsHOG(allstring1, allstring2, r2, datapath):
    # protein1 protein2 neighborhood fusion cooccurence coexpression experimental database textmining combined_score
    refs = ['neighborhood', 'fusion', 'cooccurence', 'coexpression', 'experimental', 'database', 'textmining', 'combined_score']
    final = {}
    with open(datapath, 'r') as stringdata:
        for combos in itertools.product(allstring1, allstring2):
            id1, id2 = combos
            # for id1, id2 in itertools.combinations(allstring1, allstring2):
            if r2.get(''.join(sorted([id1, id2]))) is not None:
                line = r2.get(''.join(sorted([id1, id2])))
                stringdata.seek(int(line), 0)
                words = stringdata.readline()
                # words = linecache.getline(datapath, int(line)).split()
                # words = stringdata.readline(int(line)).split()
                print(words)
                row_dict = {refs[i]: int(entry) for i, entry in enumerate(words[2:])}
                final[''.join(sorted([id1, id2]))] = row_dict
        return final


# map HOGS to OMA ID

# map OMIDS to UNIPROT

# MAP UNIPROT to string

# Find STRING vs STRING in files
