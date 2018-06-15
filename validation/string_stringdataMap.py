import pandas
import dask
import h5py
import pyoma
from utils import hashutils
import itertools

def connect2IDmap():
	r1 = redis.StrictRedis(host='localhost', port=6379, db=1)
	return r1
def connect2Stringmap():
	r2 = redis.StrictRedis(host='localhost', port=6379, db=2)
	return r2
def fam2stringID(dbobj,fam, r1):
    members = dbobj.iter_members_of_hog_id( hashutils.fam2hogID(fam)  )
    allstring = [ r.get(m.canonicalid) for m in members ]
    return allstring
def HOGvsHOG( allstring1,allstring2, r2 , datatpath):
	#protein1 protein2 neighborhood fusion cooccurence coexpression experimental database textmining combined_score
	refs = ['neighborhood', 'fusion', 'cooccurence',
	'coexpression', 'experimental', 'database', 'textmining', 'combined_score']
	final = {}
	with open( datapath , 'r')as stringdata:
		for id1, id2 in itertools.combinations(allstring1, allstring2):
			if r2( ''.join(sorted([id1, id2]) ) is not None:
				line = r2.get(''.join(sorted([id1, id2]) )
				words = stringdata.readline(line).split()
				rowdict = { refs[i]:int(entry) for i,entry in enumerate(words[2:]) }
				final[''.join(sorted([id1, id2])]= rowdict
		return final

#map HOGS to OMA ID

#map OMIDS to UNIPROT

#MAP UNIPROT to string


#Find STRING vs STRING in files
