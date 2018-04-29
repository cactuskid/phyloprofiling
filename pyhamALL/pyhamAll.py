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

from datasketch import MinHashLSH , MinHashLSHForest

import h5sparse
from scipy.sparse import vstack
import pandas as pd
from distributed import Client, progress, LocalCluster, Lock

from multiprocessing.pool import ThreadPool
if __name__ == '__main__':

	dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']
	# TODO change chucksize for full run
	chunksize = 1000
	parallel = True
	OVERWRITE = True
	saveLSH = True
	if OVERWRITE == True:
			writemode = 'w'
	else:
			writemode = 'a'

	#open up OMA
	h5OMA = open_file(config.omadir + 'OmaServer.h5', mode="r") 
	runName = time.clock()

	#setup db objects
	dbObj = db.Database(h5OMA)
	omaIdObj = db.OmaIdMapper(dbObj)
	with  h5py.File(config.datadir+ 'hashes.h5',writemode, libver = 'latest') as h5hashes:
		with h5sparse.File(config.datadir + "matrix.h5",writemode) as h5matrix:
			#with open(config.datadir + 'newlsh.pkl' , 'wb') as lshout:
				#with open(config.datadir + 'newlshforest.pkl' , 'wb') as forestout:
					#with open('errors.txt' , 'w')as errfile:

			dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']
			for dataset_name in dataset_names:
					if dataset_name not in list(h5hashes.keys()):
						dataset = h5hashes.create_dataset(dataset_name, (chunksize,0), maxshape=(None, None), dtype = 'int32')
			dsets = {}
			matrixdsets = {}
			for dataset_name in list(h5hashes.keys()):
				dsets[dataset_name] = h5hashes[dataset_name]
			for dataset_name in list(h5matrix.h5f.keys()):
				if dataset_name == 'fam':
					matrixdsets[dataset_name] = h5matrix.h5f[dataset_name]
				else:
					matrixdsets[dataset_name] = h5matrix[dataset_name]

			# corrects species tree and replacement dictionary for orthoXML files
			dic, tree = format_files.create_species_tree(h5OMA, omaIdObj)
			#set up tree mapping dictionary
			taxaIndex,reverse  = profileGen.generateTaxaIndex(tree)
			#initialize matrix h5 file

			start = time.clock()
			lsh = MinHashLSH(threshold=0.7, num_perm=128)
			forest = MinHashLSHForest(num_perm=128)
			startfam = 0
			count = 0 

			if parallel == True:
				
				
				from dask import dataframe as ddf
				import dask

				#dask.set_options(pool=ThreadPool( mp.cpu_count() ))
				#dask.set_options(get=dask.threaded.get)
				print('init cluster')
				cluster = LocalCluster(n_workers=int(mp.cpu_count()/2)  ) 

				print('init client')
				c = Client(cluster)
				print(chunksize)

				print('config functions')
				HAMPIPELINE = functools.partial( pyhamPipeline.runpyham , species_tree=tree )
				HASHPIPEline = functools.partial( profileGen.DFTree2Hashes  )
				ROWPIPELINE = functools.partial( profileGen.Tree2mat , taxaIndex = taxaIndex)
				
				print('start!')
				fams = {}
				hashdict = {}
				errors = []
				try:
					for i,fam in enumerate( pyhamPipeline.yieldFamilies(h5OMA,startfam)):
						fams[fam] = { 'ortho':pyhamPipeline.readortho( fam ,   dbObj= dbObj , species_tree=tree , replacement_dic= dic)}
						if len(fams)>chunksize:
							print(time.clock()-start)
							pddf = pd.DataFrame.from_dict(fams, orient= 'index' )	
							pddf['fams'] = pddf.index
							df = ddf.from_pandas(pddf , chunksize= 300 )
							df['tree'] = df[['fams','ortho']].apply( HAMPIPELINE , axis =1 ,  meta=pd.Series(dtype=object ) ).compute()
							hashes = df[['fams','tree']].apply( HASHPIPEline , axis =1 ,  meta=pd.Series(dtype=object) ).compute().to_dict()
							#df['rows'] = df['tree'].apply( ROWPIPELINE, meta=pd.Series(dtype=object ) ).compute()
							print(time.clock()-start)
							"""print('lsh')
																									for fam in hashes:
																											for famhashname in hashes[fam]['dict']:
																												try:
																													lsh.insert(famhashname , hashes[fam]['dict'][famhashname])
																													forest.add(famhashname ,hashes[fam]['dict'][famhashname])
																												except:
																													print('forest/lsh')
																									print(time.clock()-start)
																			
																									if saveLSH == True and time.clock() -start >2000:
																										
																										
																										print('saving')
																										pickle.dump(lsh, lshout, -1)
																										pickle.dump(forest, forestout, -1)
																			
																										start = time.clock()
																			"""
							for fam in hashes:
								try:
									hashvals = hashes[fam]['hashes']
									for i,event in enumerate(hashvals):
										dataset = dataset_names[i]
										if len(dsets[dataset])< fam+10:
											dsets[dataset].resize( (len(dsets[dataset]) +chunksize , len(hashvals[event]) ) )
										dsets[dataset][fam,:] = hashvals[event]
								except:
									print('error')
									print(fam)
									errfile.write(str(entry) + '\n')
							fams = {}
							count += chunksize
							print(count)
							print(time.clock()-start)
							print('donechunk')

				except:
					print('fam error')
				"""
														if saveLSH == True:
															pickle.dump(lsh, lshout, -1)
															pickle.dump(forest, forestout, -1)
														
													if saveLSH == True:
														forest.index()
														pickle.dump(lsh, lshout, -1)
														pickle.dump(forest, forestout, -1)
				"""

						

				print('DONE')


							
						 


	if parallel == False:


		
		

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
