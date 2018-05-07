import config
import pyhamPipeline
import profileGen
import h5sparse
from sklearn import covariance
import pyoma
import numpy as np 
import pickle
import tables
from pyoma.browser import db 
import h5py
import glob
import multiprocessing as mp
import format_files
import functools
import time
from distributed import Client, progress, LocalCluster, Lock
import pandas as pd
import dask.dataframe as ddf


if __name__ == '__main__':

	min_members = 10
	genrandqueries = 200

	clusters = []
	used_queries = set([])
	chunksize = mp.cpu_count()*4



	#open up OMA
	h5OMA = tables.open_file(config.omadir + 'OmaServer.h5', mode="r")
	#setup db objects
	dbObj = db.Database(h5OMA)
	omaIdObj = db.OmaIdMapper(dbObj)
	print(config.datadir)
	print('globbing')
	print(glob.glob(config.datadir + '*.h5'))
	dic, tree = format_files.create_species_tree(h5OMA, omaIdObj)
	
	#set up tree mapping dictionary
	taxaIndex,reverse  = profileGen.generateTaxaIndex(tree)
	
	HAMPIPELINE = functools.partial( pyhamPipeline.runpyham , species_tree=tree )
	HASHPIPEline = functools.partial( profileGen.DFTree2Hashes  )
	ROWPIPELINE = functools.partial( profileGen.Tree2mat , taxaIndex = taxaIndex)
	
	queries = list(np.random.randint(low=1, high=len(h5OMA.root.OrthoXML.Index) , size=genrandqueries))
	print(queries[0:100])

	
	dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']
	

	clusters = {}
	for cutoff in [.5, .6,.7,.8,.9]:
			clusters[cutoff] = []
	
	fams = {}
	#grab rand clusters
	for fam in queries:
		#todo limit min hog size with pyoma
		if fam not in used_queries: 
			fams[fam] = { 'ortho':pyhamPipeline.readortho( fam ,   dbObj= dbObj , species_tree=tree , replacement_dic= dic)}
		

		start = time.clock()
		print('init cluster')
		cluster = LocalCluster(n_workers=int(mp.cpu_count()/2))
		print('init client')
		c = Client(cluster)
		
		print(time.clock()-start)
		
		qpddf = pd.DataFrame.from_dict(fams, orient= 'index' )	
		qpddf['fams'] = qpddf.index
		qdf = ddf.from_pandas(qpddf , npartitions = mp.cpu_count() )
		qdf['tree'] = qdf[['fams','ortho']].apply( HAMPIPELINE , axis =1 ,  meta=pd.Series(dtype=object ) ).compute()
		qdf['hashes'] = qdf[['fams','tree']].apply( HASHPIPEline , axis =1 ,  meta=pd.Series(dtype=object) ).compute()
		hashes = qdf['hashes'].compute().to_dict()

		print('loading lsh objs')
		with open(config.datadir + 'newlsh.pkl' , 'rb') as lshout:
			lsh = pickle.load(lshout)

		with open(config.datadir + 'newlshforest.pkl' , 'rb') as lshout:
			forest = pickle.load(lshout)


		print('done')
		forestRes={}
		lshRes={}
		clustersizes = []
		for qfam in hashes:
			if hashes[qfam] is not None:
				
				results = forest.query(list(hashes[qfam]['dict'].values())[0], 200)
				print(results)

				forestRes[qfam]= set([res.split('-')[0] for res in results])
				lshresults = lsh.query(list(hashes[qfam]['dict'].values())[0])
				lshRes[qfam] = set([res.split('-')[0] for res in results])
				
				if lshRes[qfam]>0:
					print(len(lshRes[qfam]))
				
				clustersizes.append(len(lshRes[qfam]))

		with open('clustersizes.pkl', 'wb')as clusterout:
			pickle.dump(clustersizes, clusterout, -1)


"""				
		rfams = len(list(forestRes.values()))
		print(len(rfams))
		rfamdict = {}
		for fam in rfams:
			rfamdict[fam] = { 'ortho':pyhamPipeline.readortho( fam ,   dbObj= dbObj , species_tree=tree , replacement_dic= dic)}

		rpddf = pd.DataFrame.from_dict(fams, orient= 'index' )	
		rpddf['fams'] = rpddf.index
		rdf = ddf.from_pandas(rpddf , npartitions = mp.cpu_count() )
		rdf['tree'] = rdf[['fams','ortho']].apply( HAMPIPELINE , axis =1 ,  meta=pd.Series(dtype=object ) ).compute()
		rdf['hashes'] = rdf[['fams','tree']].apply( HASHPIPEline , axis =1 ,  meta=pd.Series(dtype=object) ).compute()
		rdf['rows'] = rdf[['fams','tree']].apply( ROWPIPELINE , axis =1 ,  meta=pd.Series(dtype=object) ).compute()
		hashes = rdf['hashes'].compute().to_dict()
		rows = rdf['rows'].compute().to_dict()
		





	
			for cutoff in [.5, .6,.7,.8,.9]:
				# add the new cluster in the list


			for f in result:
				if f not in used_queries:
						used_queries.update(f)
				except:
					print('error fam :' + fam)

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
"""