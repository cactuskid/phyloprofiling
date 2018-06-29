import pyham
from pyoma.browser import db 
import numpy as np
from tables import *
import re
from ete3 import Tree
import pickle
import tempfile
import functools
import config
import h5py
import time
import pyhamPipeline
import profileGen
import format_files
import multiprocessing as mp

import functions
import datetime

from datasketch import MinHashLSH , MinHashLSHForest

import h5sparse
from scipy.sparse import vstack
import pandas as pd
from distributed import Client, progress, LocalCluster, Lock

from dask import dataframe as ddf
import dask

from multiprocessing.pool import ThreadPool
if __name__ == '__main__':
	dt = 
	datestr = "{%B%d%Y_%H%M}".format(datetime.now())
	print(datestr)

	dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']
	# TODO change chucksize for full run
	chunksize = 1000
	parallel = True
	testmode = True

	#open up OMA
	h5OMA = open_file(config.omadir + 'OmaServer.h5', mode="r") 
	dbObj = db.Database(h5OMA)
	omaIdObj = db.OmaIdMapper(dbObj)

		
	

	treemap_fam = pyhamPipeline.get_hamTree(1, dbObj, tree, dic)
	matrixRow = profileGen.Tree2mat(treemap_fam, taxaIndex)
	
	for dataset_name in ['fam', 'hogmatrix']:
		if dataset_name not in list(h5matrix.h5f.keys()):
			if dataset_name != 'fam':
				dataset = h5matrix.create_dataset(dataset_name, data=matrixRow , chunks=(100000,), maxshape=(None,))
			else:
				dataset = h5matrix.h5f.create_dataset(dataset_name, (chunksize,0), maxshape=(None, None), dtype = 'int32')


	
	#temporary, 
	rows = []
	errors = []

	#load Fam
	for i,fam in enumerate(pyhamPipeline.yieldFamilies(h5OMA, startfam)):
		
		#generate treemap profile
		try:
			treemap_fam = pyhamPipeline.get_hamTree(fam, dbObj, tree, dic)
			
			# generate matrix of hash
			hashesDic = profileGen.Tree2Hashes(treemap_fam, fam, lsh)
			matrixRow = profileGen.Tree2mat(treemap_fam, taxaIndex, verbose = False )
			rows.append(matrixRow)

			if i == 0:
				for dataset in dsets:
					if dataset == 'fam':
						dsets[dataset].resize((startfam+ chunksize+1, 1 ))
					else: 
						dsets[dataset].resize((startfam + chunksize+1, np.asarray(hashesDic[dataset]).shape[0] ))
						
			if i % chunksize == 0 and i != 0:
				print(time.clock()-start)
				with open(config.datadir + 'lsh.pkl' , 'wb') as lshout:
					pickle.dump(lsh, lshout, -1)
				
				print(fam)
				print(i)

				if time.clock()-start > 1000:
					with open(config.datadir + 'lsh.pkl' , 'wb') as lshout:
						pickle.dump(lsh, lshout, -1)
					start = time.clock()


				for dataset in dsets:
					if dataset == 'fam':
						dsets[dataset].resize((startfam + i+chunksize+1, 1 ))
					else: 
						dsets[dataset].resize((startfam + i+chunksize, np.asarray(hashesDic[dataset]).shape[0] ))
				
				matrixdsets['hogmatrix'].append(vstack(rows))
				rows= []
			for dataset in dsets:
				if dataset == 'fam':
					dsets[dataset][i,:] = fam
				else:
					dsets[dataset][i,:] = hashesDic[dataset]
			
		except:
			errors.append(fam)
			with open(config.datadir + 'errors.pkl' , 'wb') as errout:
				pickle.dump(errors, errout, -1)
			
	with open(config.datadir + 'lsh.pkl' , 'wb') as lshout:
		pickle.dump(lsh, lshout, -1)
	

	with open(config.datadir + 'errors.pkl' , 'wb') as errout:
		pickle.dump(errors, errout, -1)
	
	print('DONE!')

	h5hashes.close()
	h5OMA.close()
