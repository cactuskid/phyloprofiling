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
from datasketch import MinHashLSH
import h5sparse
from scipy.sparse import vstack
import pandas as pd


dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']
# TODO change chucksize for full run
chunksize = 1000
parallel = True
OVERWRITE = False

if OVERWRITE == True:
		writemode = 'w'
	else:
		writemode = 'a'

#open up OMA
h5OMA = open_file(config.omadir + 'OmaServer.h5', mode="r") 

	

h5hashes = h5py.File(config.datadir+ 'hashes.h5',writemode)
h5matrix =  h5sparse.File(config.datadir + "matrix.h5",writemode)

dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']




startfam = 0
if OVERWRITE ==False:
	startfam = np.amax(h5hases['fam'])

#setup db objects
dbObj = db.Database(h5OMA)
omaIdObj = db.OmaIdMapper(dbObj)


#load Fams
if parallel == True:
	
	
	from dask import dataframe as ddf
	from dask.distributed import Client
	c = Client()
	l = distributed.Lock(name='OMAlock', client=c)

	df = ddf.read_hdf(config.omadir + 'OmaServer.h5', '/root/OrthoXML/Index')
	
	for line in pyhamPipeline.yieldFamilies(h5file,startfam)
		fams.append(fam)
		if len(fams)>chunksize:
			df = ddf.from_array()

			fams = []

	print(len(df))
	
	"""pint(df)
			
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
			"""



if parallel == False:


	for dataset_name in dataset_names:
		if dataset_name not in list(h5hashes.keys()):
			dataset = h5hashes.create_dataset(dataset_name, (chunksize,0), maxshape=(None, None), dtype = 'int32')


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

	
	#temporary, 
	rows = []
	errors = []


	start = time.clock()
	lsh = MinHashLSH()
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
						dsets[dataset].resize((chunksize+1, 1 ))
					else: 
						dsets[dataset].resize((chunksize+1, np.asarray(hashesDic[dataset]).shape[0] ))
						
			if i % chunksize == 0 and i != 0:
				print(time.clock()-start)
				with open(config.datadir + 'lsh.pkl' , 'wb') as lshout:
					pickle.dump(lsh, lshout, -1)
				
				print(fam)
				print(i)

				for dataset in dsets:
					if dataset == 'fam':
						dsets[dataset].resize((i+chunksize+1, 1 ))
					else: 
						dsets[dataset].resize((i+chunksize, np.asarray(hashesDic[dataset]).shape[0] ))
				
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
