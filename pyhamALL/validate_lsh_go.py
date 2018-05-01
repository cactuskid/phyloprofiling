import tables
import numpy as np
import pickle

from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import TermCounts 
from pyoma.browser import db

import config

import validation.phyloValidationGoTerm
from profileGen import get_HOGhash, jaccard_rank

from datasketch import MinHashLSH , MinHashLSHForest
import h5py


# load all the files (oma database, go dag file, association gaf file)
# set oma database
# contains genes and related GO terms
h5file = tables.open_file(config.omadir + 'OmaServer.h5', mode='r')

dbObj = db.Database(h5file)
omaIdObj = db.OmaIdMapper(dbObj)
#load gene ontology DAG
#use guaranteed acyclic basic GO
# TODO provide correct path for those
go = obo_parser.GODag('/home/laurent/Documents/phylo/phyloprofiling_working/1stdraft/data/go.obo')
associations = read_gaf('/home/laurent/Documents/phylo/phyloprofiling_working/1stdraft/data/gene_association.tair.gz')

# Get the counts of each GO term.
termcounts = TermCounts(go, associations)

goTermAnalysis = phyloValidationGoTerm.SemanticSimilarityAnalysis(go, h5file, termcounts)

#### LSH cluster

queries = []
#clustering
# list of list
clusters = []
# list of used fam
used_queries = set([])

lsh = pickle.lead(open(config.datadir + 'lsh.pkl' , 'rb'))

import multiprocessing as mp

pool = mp.Pool()


def jaccardmp( h1 , h2):
	i,h1 = h1
	j,h2 = h2

	return( i,j, h1.jaccard(h2) )


with  h5py.File(config.datadir+ 'hashes.h5', 'r') as h5hashes:

	genrandqueries = 100
	#np.random.seed(1)
	queries = list(np.random.randint(low=1, high=len(h5OMA.root.OrthoXML.Index) , size=genrandqueries))
	#filter . if has go terms for gene function	
	results = {}
	for combo in itertools.combinations(['duplication', 'gain', 'loss', 'presence'])
		for fam in fam_list_queries:
			
			query = profileGen.gethoghash(fam, combo)	
			reuslts = (lsh.query(query)
			
			#filter . if has go terms for gene function

			hashes = [ profileGen.gethoghash(hog, combo)  for hog in results +[fam]]
			go = [ getgoterms(hog) for hog in results+[fam] ]
			results[fam] = { 'fam' : results +[fam] , 'hash': hashes , 'go': go  }

			godist = np.zeros( (len(go), len(go)))
			for i,go1 in enumerate(go):
				for j,go2 in enumerate(go):
					godist[i,j] = dist(go1,go2)

			jdist = np.zeros( (len(go), len(go)))
			hashcompare = [ (h1,h2) for h1,h2 in itertools.combinations(enumerate(hashes), 2) ]
			results = pool.map_async(hashcompare , jaccardmp ).get()
			for r in results:
				i,j,dist = r
				jdist[i,j] = dist



	#validate stuff
	for fam in fam_list_quries:
		for result in results[fam]:
			results_go[fam][result] = goTermAnalysis.semantic_similarity_score(fam, result)

	# plot if same
