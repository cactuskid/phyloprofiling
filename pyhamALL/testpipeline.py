import config
import pyhamPipeline
import profileGen
import h5sparse
from sklearn import Covariance
import pyoma
import numpy as np 


min_members = 10
genrandqueries = 100

clusters = []
used_queries = set([])

queries = list(np.randint(low=1, high=10000000000, size=genrandqueries))

#load lsh
with open(config.datadir + 'newlsh.pkl' , 'rb') as lshout:
	lsh = pickle.load(lshout)

#open up OMA
h5OMA = open_file(config.omadir + 'OmaServer.h5', mode="r") 

#setup db objects
dbObj = db.Database(h5OMA)
omaIdObj = db.OmaIdMapper(dbObj)
	

h5hashes = h5py.File(config.datadir+ 'hashes.h5',writemode)
h5matrix =  h5sparse.File(config.datadir + "matrix.h5",writemode)
dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']

def get_result(query, raw_result):
	filtered_results = [query] + [x for x in raw_result if x.split('-')[1] == query.split('-')[1]]
	return filtered_results

clusters = {}
for cutoff in [.5, .6,.7,.8,.9]:
		clusters[cutoff] = []



#grab rand clusters
with open('clustersize.txt', 'a')as clustersizeout:
	for fam in queries:
		#todo limit min hog size with pyoma
		
		if fam not in used_queries:       
			queryhash = get_HOGhash(fam , h5hashes) 
			result = get_result(query, lsh.query(queryhash))
			rhashes= {}
			for rfam in result:
				rhashes[rfam]= getHOGhash(rfam)
			fams, scores = jaccard_rank(queryhash , rhashes )
			for cutoff in [.5, .6,.7,.8,.9]:
				# add the new cluster in the list
				clusters[cutoff].append(jaccard_cutoff(fams,scores,cutoff))
				clustersizeout.write( str(cutoff) len())
				# exlude the results from the fams         
				for f in result:
					if f not in used_queries:
						used_queries.update(f)


with open(config.datadir +  'rundata.txt' , 'a')as rundata:
	times = []
	sizes = []
	for cutoff in clusters:
		for results in clusters[cutoff]:
			start = time.clock()
			rows = []
			#todo generate hogmat on the fly

			#run single core
			model = GraphLassoCV(verbose=True, njobs=1)
			model.fit(hogmat)
			cov_ = model.covariance_
			times = time.clock() - start
			size = len(c)
			rundata.write(str(time)  +','+ str(size)+'\n')
