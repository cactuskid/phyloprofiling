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
from datetime import datetime

from datasketch import MinHashLSH , MinHashLSHForest

import h5sparse
from scipy.sparse import vstack
import pandas as pd

from distributed import Client, progress, LocalCluster, Lock, Scheduler, Adaptive, Variable

from dask import dataframe as ddf
from dask import delayed


from multiprocessing.pool import ThreadPool
if __name__ == '__main__':
	
	
	#open up OMA

	with open_file(config.omadir + 'OmaServer.h5', mode="r") as h5OMA:
		dbObj = db.Database(h5OMA)
		omaIdObj = db.OmaIdMapper(dbObj)
		dic, tree = format_files.create_species_tree(h5OMA, omaIdObj)
		taxaIndex,reverse  = profileGen.generateTaxaIndex(tree)
		HAMPIPELINE = functools.partial( pyhamPipeline.runpyham , species_tree=tree )
		HASHPIPEline = functools.partial( profileGen.DFTree2Hashes  )
		ROWPIPELINE = functools.partial( profileGen.Tree2mat , taxaIndex = taxaIndex)
		

	def genDFs( h5OMA, size = 100):
		fams = {}
		with open_file(config.omadir + 'OmaServer.h5', mode="r") as h5OMA:
			dbObj = db.Database(h5OMA)
			omaIdObj = db.OmaIdMapper(dbObj)
			dic, tree = format_files.create_species_tree(h5OMA, omaIdObj)
			taxaIndex,reverse  = profileGen.generateTaxaIndex(tree)
			READORTHO = functools.partial( pyhamPipeline.readortho , dbObj= dbObj , species_tree=tree , replacement_dic= dic )
			for i,row in enumerate( h5OMA.root.OrthoXML.Index):
				fam = row[0]
				fams[fam] = { 'ortho':READORTHO( fam ) }
				if len(fams)>size: 
					pddf = pd.DataFrame.from_dict(fams, orient= 'index' )
					pddf['Fam'] = pddf.index
					fams = {}
					yield pddf

	def worker(i,q, retq,l):
		print('worker init ' + str(i))
		while True:
			df = q.get()
			if df is not None:
				df['tree'] = df[['Fam','ortho']].apply( HAMPIPELINE , axis =1  )
				df['hash'] = df[['Fam','tree']].apply( HASHPIPEline , axis =1 )
				#df['rows'] = df['tree'].apply( ROWPIPELINE )
				retq.put(df)
			else:
				print('DONE WORKER'+ str(i))
				break


	def saver( i,q,retq,l ):
		
		dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']

		printstart = time.clock()
		savestart = time.clock()
		globaltime = time.clock()
		
		
		chunksize = 100
		count = 0 
		dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']
		threshold=0.7
		lsh = MinHashLSH(threshold=threshold, num_perm=128)
		forest = MinHashLSHForest(num_perm=128)

		dt = datetime
		datestr = "{:%B_%d_%Y_%H_%M}".format(dt.now())

		with open('errors.txt' , 'w')as errfile:
			
			with  h5py.File(config.datadir+ datestr + 'hashes.h5','w', libver = 'latest') as h5hashes:
				dsets = {}
				for dataset_name in dataset_names:
					if dataset_name not in list(h5hashes.keys()):
						dataset = h5hashes.create_dataset(dataset_name, (chunksize,0), maxshape=(None, None), dtype = 'int32')
					dsets[dataset_name] = h5hashes[dataset_name]
				
				with h5sparse.File(config.datadir + datestr + "matrix.h5",'w'  ) as h5matrix:
					with open(config.datadir + datestr + '_'+str(threshold)+'_' +'newlsh.pkl' , 'wb') as lshout:
						
						with open(config.datadir +datestr+ 'newlshforest.pkl' , 'wb') as forestout:	
							print('saver init ' + str(i))
							while True:
								thisdf = retq.get()
								print(time.clock()- globaltime)
								print(count)
								if thisdf is not None:
									
									hashes = thisdf['hash'].to_dict()

									#rows = vstack(thisdf['rows'])
									
									for fam in hashes:
										

										if hashes[fam] is not None:

											for famhashname in hashes[fam]['dict']:
												lsh.insert(famhashname , hashes[fam]['dict'][famhashname])
												forest.add(famhashname ,hashes[fam]['dict'][famhashname])
											hashvals = hashes[fam]['hashes']
											
											for i,event in enumerate(hashvals):
												dataset = dataset_names[i]
												if len(dsets[dataset])< fam+10:
													dsets[dataset].resize( (fam +chunksize , len(hashvals[event])  ) )
												dsets[dataset][fam,:] = hashvals[event]
											h5hashes.flush()
										else:
											print('error')
											print(fam)
											errfile.write(str(fam) + '\n')


										
									if time.clock() - printstart > 60:
										print(thisdf['Fam'].max())	
										print( time.clock()- globaltime)
										printstart = time.clock()
									
									if time.clock() -savestart > 2000:
										print(thisdf['Fam'].max())
										print('saving')
										
										pickle.dump(lsh, lshout, -1)
										pickle.dump(forest, forestout, -1)
										
										savestart = time.clock() 
										print( time.clock()- globaltime)
									count+= len(thisdf)
								else:
									pickle.dump(lsh, lshout, -1)
									pickle.dump(forest.index(), forestout, -1)
									
									print('DONE UPDATER' + str(i))
									break

	def mp_with_timeout(nworkers, nupdaters, datagenerator , workerfunction, updatefunction ):
		wprocesses ={}
		uprocesses ={}    
		l = mp.Lock()
		cores = mp.cpu_count()
		q = mp.Queue( maxsize = cores*10 )
		retq = mp.Queue( maxsize = cores*10 )
		

		print('start workers')
		for i in range(nworkers):
			t = mp.Process(target=workerfunction, args=(i,q,retq,l)  ) 
			t.daemon = True
			t.start()
			wprocesses[i] = t
		print('start savers')
		for i in range(nupdaters):
			t = mp.Process(target=updatefunction, args=(i,q,retq,l) ) 
			t.daemon = True
			t.start()
			uprocesses[i]= t

		count =0
		start = time.time()
		while True:
			time.sleep(.1)
			try:
				data = next(datagenerator)
				q.put(data)
				count += 1
			except StopIteration:	
				print( 'stop iteration')
				for p in range(nworkers):
					q.put(None)
				for p in range(nupdaters):
					retq.put(None)
				break
		for p in wprocesses:
			wprocesses[p].terminate()
		gc.collect()
		for p in uprocesses:
			uprocesses[p].terminate()
		gc.collect()
		print( 'DONE!!!!!')

	mp_with_timeout(int(mp.cpu_count()/2), 1 , genDFs(h5OMA,100) , worker, saver )