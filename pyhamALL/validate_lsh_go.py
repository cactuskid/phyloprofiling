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

with  h5py.File(config.datadir+ 'hashes.h5', 'r') as h5hashes:

	fam_list_queries = [23, 4345, 2134, 25]

	fam_hashdict = get_hashdict(fams, h5hashes)
	
	results = {}

	for fam in fam_list_queries:
		results[fam](lsh.query(fam_hashdict[fam])
	
	# call jaccard_rank ??

	#validate stuff
	for fam in fam_list_quries:
		for result in results[fam]:
			results_go[fam][result] = goTermAnalysis.semantic_similarity_score(fam, result)

	# plot if same
