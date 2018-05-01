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
	
	dt = datetime
	
	datestr = "{:%B_%d_%Y_%H_%M}".format(dt.now())
	print(datestr)

	dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']
	# TODO change chucksize for full run
	chunksize = 10
	parallel = True
	testmode = True

	#open up OMA
	h5OMA = open_file(config.omadir + 'OmaServer.h5', mode="r") 
	dbObj = db.Database(h5OMA)
	omaIdObj = db.OmaIdMapper(dbObj)
	
	start = time.clock()
	lsh = MinHashLSH(threshold=0.7, num_perm=128)
	forest = MinHashLSHForest(num_perm=128)
	startfam = 0
	count = 0 

	dic, tree = format_files.create_species_tree(h5OMA, omaIdObj)
	taxaIndex,reverse  = profileGen.generateTaxaIndex(tree)
	dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']


	print('init cluster')

	cluster = LocalCluster(n_workers=1 ) 
	scheduler = Scheduler()
	#adapative_cluster = Adaptive(scheduler, cluster)
	scheduler.start()

	print('init client')
	c = Client(cluster)

	print('config functions')
	READORTHO = functools.partial( pyhamPipeline.readortho , dbObj= dbObj , species_tree=tree , replacement_dic= dic )
	HAMPIPELINE = functools.partial( pyhamPipeline.runpyham , species_tree=tree )
	HASHPIPEline = functools.partial( profileGen.DFTree2Hashes  )
	ROWPIPELINE = functools.partial( profileGen.Tree2mat , taxaIndex = taxaIndex)

	row0= ROWPIPELINE(  HAMPIPELINE((1,READORTHO(1))) )
	print(row0)
	
	@delayed(pure = True)
	def retDFs(start,stop, h5OMA):
		fams = {}
		for i,fam in enumerate(range(start,stop)):
			fams[fam] = { 'ortho':READORTHO( fam ) } 
		pddf = pd.DataFrame.from_dict(fams, orient= 'index' )	
		pddf['fams'] = pddf.index
		fams = {}
		return pddf


	with  h5py.File(config.datadir+ datestr + 'hashes.h5','w', libver = 'latest') as h5hashes:
		with h5sparse.File(config.datadir + datestr + "matrix.h5",'w'  ) as h5matrix:
			h5matrix.create_dataset('rows', data=row0)
			with open(config.datadir + datestr + 'newlsh.pkl' , 'wb') as lshout:
				with open(config.datadir +datestr+ 'newlshforest.pkl' , 'wb') as forestout:
					with open('errors.txt' , 'w')as errfile:
						#init hash h5
						dsets = {}
						for dataset_name in dataset_names:
								if dataset_name not in list(h5hashes.keys()):
									dataset = h5hashes.create_dataset(dataset_name, (chunksize,0), maxshape=(None, None), dtype = 'int32')
								dsets[dataset_name] = h5hashes[dataset_name]
						
						dataframes = [retDFs(start*chunksize,(start+1)*chunksize) for start in range(int(len(h5OMA.root.OrthoXML.Index)/chunksize) +1) ]
						df = ddf.from_delayed(dataframes , meta = {'fams':int, 'ortho':str  } )

						df['tree'] = df[['fams','ortho']].apply( HAMPIPELINE , axis =1 ,  meta=pd.Series(dtype=object ) )
						df['hash'] = df[['fams','tree']].apply( HASHPIPEline , axis =1 ,  meta=pd.Series(dtype=object) )
						df['rows'] = df['tree'].apply( ROWPIPELINE, meta=pd.Series(dtype=object ) )
						
						dfs = df.to_delayed()
						
						@delayed(pure = False)
						def save( thisdf , lsh , forest , dsets , h5matrix , lshout , forestout):
							hashes = thisdf['hash'].compute().to_dict()
							rows = vstack(thisdf['rows'].compute())
							savestart = Variable('savestart')
							printstart = Variable('printstart')
							globaltime = Variable('globalstart')

							with Lock('outputlock'):
								h5matrix.append(rows)
								for fam in hashes:
									if hashes[fam] is not None:
										lsh.insert(famhashname , hashes[fam]['dict'][famhashname])
										forest.add(famhashname ,hashes[fam]['dict'][famhashname])
										
										hashvals = hashes[fam]['hashes']
										for i,event in enumerate(hashvals):
											dataset = dataset_names[i]
											if len(dsets[dataset])< fam+10:
												dsets[dataset].resize( (len(dsets[dataset]) +chunksize , len(hashvals[event]) ) )
											dsets[dataset][fam,:] = hashvals[event]								
										else:
											print('error')
											print(fam)
											errfile.write(str(entry) + '\n')
							
							if time.clock() - printstart.get() > 60:
								with Lock('outputlock'):
									print(thisdf['fams'].max())	
									print( time.clock()- globaltime.get())

								printstart.set(time.clock())

							if time.clock() -savestart.get() >2000:
								
								with Lock('outputlock'):
									print(thisdf['fams'].max())
									print('saving')
									pickle.dump(lsh, lshout, -1)
									pickle.dump(forest, forestout, -1)
									start.set( time.clock() )
									print( time.clock()- globaltime.get())
						
						start = Variable('savestart') 
						printstart = Variable('printstart')
						globaltime = Variable('globalstart')
						printstart.set( time.clock() )
						start.set( time.clock() )
						globaltime.set( time.clock() )
						
						SAVE = functools.partial(save , lsh=lsh , forest=forest , dsets=dsets , h5matrix=h5matrix , lshout=lshout , forestout=forestout )
						print('start!')
						writes = [SAVE(df) for df in dfs ]
						ddf.compute(*writes)
						cluster.close()
						print('DONE')
