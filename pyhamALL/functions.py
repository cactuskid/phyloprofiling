

import csv
import numpy as np
import glob
import functools
import pickle
import random
import itertools
import re, string
from subprocess import Popen
import os
import multiprocessing as mp
import pandas as pd
import shlex, subprocess
import config
import dask.dataframe as dd
import dask.array as da
import dask
from dask.delayed import delayed
import h5py
import gc

dask.set_options(get=dask.threaded.get)

########################hdf5 save and load #################################################################
def hd5save(df, name , overwrite ,verbose = False , pipelines= None ):
	#dataframe columns should all be arrays or bytestrings encoding
	if overwrite == True:
		f = h5py.File(name,'w')
	else:
		f = h5py.File(name,'a')
	f.create_group('datasets')

	if pipelines is None:
		pipelines = df.columns
	
	total = list(set(pipelines+['Y']))
	for col in total:
		try:
			array = np.vstack( df[col].values )
		except:
			maxlen = 0
			for array in df[col].values:
				if len(array) > maxlen :
					maxlen=len(array)
			#pad bytestrings with spaces to maxlen

			array = np.vstack( [ np.string_(np.pad(array,((0, maxlen - len(array))) , mode='constant' , constant_values=20 )) for array in df[col].values ]  )
		
		if col not in f['datasets']:
			try:
				dset = f['datasets'].create_dataset(col, data=array ,chunks=True)
			except :
				dset = f['datasets'].create_dataset(col, data=array,  dtype="S" , chunks=True)
		else:
			dset = f['datasets'][col]
			x,y = dset.shape
			inx,iny = array.shape
			#resize dataset for new data.
			dset.resize(inx+x, y + max(0,iny-y) )
			dset[x:inx + x , : ] = array
	f.close()

def DaskArray_hd5loadDataset(files , verbose = False ):
	#load to dask arrays, all arrays passed should have the same dataset names
	datasets = {}
	for name in files:
		f = h5py.File(name,'r')
		print(list(f['datasets'].keys()) )
		for dataset in f['datasets'].keys():
			chunksize = [  max(1, int(chunk/2) ) for chunk in f['datasets'][dataset].chunks ]
			
			if verbose == True:
				print('chunking smaller than hdf5')
				print( chunksize)
				print( f['datasets'][dataset].chunks)
				#print(f['datasets'][dataset][0:10])
			if dataset not in datasets:
				datasets[dataset] = da.from_array(f['datasets'][dataset], chunks=chunksize )
				if verbose==True:
					f['datasets'][dataset][0:2]
					print(datasets[dataset][0:2].compute(get=dask.get) )
			else:
				array = da.from_array(f['datasets'][dataset], chunks=chunksize )
				datasets[dataset] = da.concatenate([array, datasets[dataset]], axis=0)
				if verbose ==True:
					print('append')
					print(datasets[dataset])
					print(datasets[dataset][0:10].compute() )
	for key in datasets:
		print(datasets[key][0:10].compute() )
	return datasets

##################################run processes, apply to dataframe partitions etc###########################################

def runOnDelayed(DF , pipeline):
	#split df into temp fastas
	dfs = DF.to_delayed()
	retdfs = [functions.delayed(pipeline) (df) for df in dfs]
	DDF =None
	for df in retdfs:
		if DDF == None:
			DDF = dd.from_pandas(df , npartitions = mp.cpu_count() )
		else:
			DDF.append(dd.from_pandas(df))
	return DDF


def openprocess(args , inputstr =None , verbose = False ):
	args = shlex.split(args)
	p = subprocess.Popen(args,  stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr= subprocess.PIPE )
	if verbose == True:
		print(inputstr.decode())
	
	if inputstr != None:
		p.stdin.write(inputstr.encode())
		
	output = p.communicate()
	if verbose == True:
		print(output)
	p.wait()
	return output[0].decode()

#######################################pipeline building functions##############################################
def retx(x):
	return x

def compose(functions):
	def compose2(f, g):
		def fOg(x):
			return f(g(x))
		return fOg
	retfunction = functools.reduce(compose2, functions, retx )
	return retfunction

###########################################################dataframe / dataset building ##########################################

def applypipeline_to_series(series, pipeline):
	print(series)
	newseries = series.map( pipeline )
	print(newseries)
	return newseries

