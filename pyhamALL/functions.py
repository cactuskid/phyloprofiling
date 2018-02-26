

import csv
import numpy as np
import scipy.signal as sig

import scipy.linalg as linalg
from scipy.fftpack import rfft, fftshift
import glob

import functools


import pickle
import random


import re, string 

from keras.models import Sequential
from keras.layers import Dense, Conv1D, Dropout
from keras.wrappers.scikit_learn import KerasClassifier
from keras.utils import np_utils
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold , train_test_split
from sklearn.preprocessing import LabelEncoder , robust_scale , normalize
from sklearn.pipeline import Pipeline
from subprocess import Popen
import os
import multiprocessing as mp
from Bio import SeqIO

import pandas as pd
import shlex, subprocess
import config
import dask.dataframe as dd
import dask.array as da
import dask
from dask.delayed import delayed
import h5py

import io


########################hdf5 save and load #################################################################
def hd5save(df, name , overwrite ,verbose = False):
    #output dataframe columns in as a matrices and store to hdf5. 

    #dataframe columns should all be arrays or bytestrings w ascii encoding

    #this will be used to save our matrix with our encoding of evolutionary events and presence over taxa.
    if overwrite == True:
        f = h5py.File(name,'w')
    else:
        f = h5py.File(name,'a')
    f.create_group('datasets')

    for col in df.columns:
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

def dfdumptohdf5(dd, name , overwrite ,verbose = False):
    #dataframe straight to hdf5
    #this will store the parts of our dataframe we don't want to represent as a matrix.
    #all strings should be in bytestr format

    #this will be used to store our pyham trees  and minhash objects.
    #name can be a pattern like 'myfile.*.hdf5' to split into multiple files

    if overwrite == True:
        dd.to_hdf(name, '/data', get=dask.multiprocessing.get) 
    else:
        #load prexisting and append
        old = dd.read_hdf(name, '/data')
        old.merge(dd, left_on='index', right_on='index', how='outer')
        old.to_hdf(name, '/data', get=dask.multiprocessing.get) 
    return name


def DaskArray_hd5loadDataset(files , verbose = False ):
    #load to dask arrays, all arrays passed should have the same dataset names
    #this will load the matrix in each HDF5 into dask arrays and concatenate them returning a big matrix with everything.

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
                    print(datasets[dataset][0:10].compute(get=dask.get) )
    return datasets



##################################run processes, apply to dataframe partitions etc###########################################

def runOnDelayed(DF , pipeline):
    #split df nd run a pipeline on each dataframe
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
    #if we need to run any external processes we can use this in the pipeline, but I doubt we'll need it
    args = shlex.split(args)
    p = subprocess.Popen(args,  stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr= subprocess.PIPE)
    if verbose == True:
        print(inputstr.decode())
    if inputstr != None:
        p.stdin.write(inputstr)
        
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


###################################### functions for phyloprofile data generation #################################


def generateTaxaIndex(newick):
    #use this to generate an index for the global taxonomic tree for all of OMA
    #pass it in hyper params

    t =ete3.Tree( newick)
    taxaIndex = {}
    for i,n in enumerate(t.traverse):
        taxaIndexReverse[i] = node.name
        taxaIndex[node.name] = i
    return taxaIndex, taxaIndexReverse





def serializeHash(minhas, hyperparams):
    #convert hash to bystr to store in HDF5

    pass 

def Ham2Hashes(ham , hyperparams):
    #encode evolutionary data in 4 hashes


    pass

def Ham2matrixRow( ham, hyperparams):
    #encode evolutionary data in 4 matrices

    pass

def serializeHam(minhas, hyperparams):
    #convert Ham tree to bystr to store in HDF5 

    pass 



def ortho2Ham(strIO,hyperparams):
    #run pyham on a single orthoxml string and output ete3 tree of evolutionary events

    pass

def orth2strIO(ortho, hyperparams):
    #generate a stringIO object to pass to pyham


    pass


def dataGenHOGs( famtable , dbObj ,hyperparams):
    #generate orthoxml strings from pyoma object
    for fam in fams:
        ortho = dbObj.get_orthoxml(fam)
        yield fam,ortho



###########################################################dataframe / dataset building ##########################################
def HogsToDF(orthoGen ,DDF = None):
    DFdict={}
    for fam,ortho in orthoGen:
        #use Natasha's namehack


        DFdict[fam] = hackedortho


    df = pd.DataFrame.from_dict(DFdict, orient = 'index')
    
    if DDF == None:
        DDF = dd.from_pandas(df , npartitions = mp.cpu_count() )
    else:
        DDF.append(dd.from_pandas(df , npartitions = mp.cpu_count() ))
    
    return DDF

