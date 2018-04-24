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

import pyhamPipeline
import profileGen
import format_files
from datasketch import MinHashLSH
import h5sparse
from scipy.sparse import vstack

dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']
# TODO change chucksize for full run
chunksize = 3
parallel = False

#load Fams
if parallel == True:
	
	
	from dask import dataframe as ddf
	from dask.distributed import Client
	c = Client()
	l = distributed.Lock(name='OMAlock', client=c)
	df = functions.famsToDF(h5file)
	#from family number to ete3
	HAMPIPELINE = functools.partial( pyhamPipeline.get_hamTree,  dbObj= dbObj , species_tree=species_tree , replacement_dic= replacement_dic, l=l)
	#from ete3 to hashes
	HASHPIPEline =  profileGen.Tree2Hashes 
	#from ete3 to matrix rows
	ROWPIPELINE = functools.partial( profileGen.Tree2mat , taxaIndex = taxa_index)

	#calculate all profiles
	df['Tree']= df['index'].map(HAMPIPELINE).compute()
	df['Hashes'] = df['Tree'].map(HASHPIPEline).compute()
	df['Rows'] = df['Tree'].map(ROWPIPELINE).compute()
	df.index.to_delayed()


if parallel == False:

	h5hashes = h5py.File(config.datadir+ 'hashes.h5','w')
	h5matrix =  h5sparse.File(config.datadir + "matrix.h5",'w')
	dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']

	for dataset_name in dataset_names:
		if dataset_name not in list(h5hashes.keys()):
			dataset = h5hashes.create_dataset(dataset_name, (chunksize,0), maxshape=(None, None), dtype = 'int32')

	#open up OMA
	h5OMA = open_file(config.omadir + 'OmaServer.h5', mode="r") 


	#setup db objects
	dbObj = db.Database(h5OMA)
	omaIdObj = db.OmaIdMapper(dbObj)

	# corrects species tree and replacement dictionary for orthoXML files
	dic, tree = format_files.create_species_tree(h5OMA, omaIdObj)
	#set up tree mapping dictionary
	taxaIndex,reverse  = profileGen.generateTaxaIndex(tree)
	#initialize matrix h5 file

	treemap_fam = pyhamPipeline.get_hamTree(1, dbObj, tree, dic)
	matrixRow = profileGen.Tree2mat(treemap_fam, taxaIndex)
	
	for dataset_name in ['fam', 'hogmatrix']:
		if dataset_name not in list(h5matrix.h5f.keys()):
			if dataset_name != 'fam':
				dataset = h5matrix.create_dataset(dataset_name, data=matrixRow , chunks=(100000,), maxshape=(None,))
			else:
				dataset = h5matrix.h5f.create_dataset(dataset_name, (chunksize,0), maxshape=(None, None), dtype = 'int32')

	dsets = {}
	matrixdsets = {}
	for dataset_name in list(h5hashes.keys()):
		dsets[dataset_name] = h5hashes[dataset_name]
		
	for dataset_name in list(h5matrix.h5f.keys()):
		if dataset_name == 'fam':
			matrixdsets[dataset_name] = h5matrix.h5f[dataset_name]
		else:
			matrixdsets[dataset_name] = h5matrix[dataset_name]

	
	chunksize = 1000
	#temporary, 
	rows = []
	lsh = MinHashLSH()
	#load Fam
	for i,fam in enumerate(pyhamPipeline.yieldFamilies(h5OMA)):
		#generate treemap profile
		treemap_fam = pyhamPipeline.get_hamTree(fam, dbObj, tree, dic)
		# generate matrix of hash
		hashesDic = profileGen.Tree2Hashes(treemap_fam, fam, lsh)
		matrixRow = profileGen.Tree2mat(treemap_fam, taxaIndex, verbose = True )
		rows.append(matrixRow)

		if i == 0:
			for dataset in dsets:
				if dataset == 'fam':
					dsets[dataset].resize((chunksize+1, 1 ))
					matrixdsets[dataset].resize((chunksize+1, 1 ))
				else: 
					dsets[dataset].resize((chunksize+1, np.asarray(hashesDic[dataset]).shape[0] ))
					
		if i % chunksize == 0 and i != 0:
			for dataset in dsets:
				if dataset == 'fam':
					dsets[dataset].resize((i+chunksize+1, 1 ))
					matrixdsets[dataset].resize((chunksize+1, 1 ))
				else: 
					dsets[dataset].resize((i+chunksize, np.asarray(hashesDic[dataset]).shape[0] ))
			print(vstack(rows))
			matrixdsets['hogmatrix'].append(vstack(rows))
			print(matrixdsets['hogmatrix'][:].toarray())
			rows= []
		for dataset in dsets:
			if dataset == 'fam':
				dsets[dataset][i,:] = fam
				matrixdsets[dataset][i,:] = fam
			else:
				dsets[dataset][i,:] = hashesDic[dataset]
		

	
	with open(config.datadir + 'lsh.pkl' , 'wb') as lshout:
		pickle.dump(lsh, lshout, -1)
	
	h5hashes.close()
	h5OMA.close()
