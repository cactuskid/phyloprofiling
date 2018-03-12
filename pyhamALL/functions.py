

import csv
import numpy as np
import scipy.signal as sig

import scipy.linalg as linalg
from scipy.fftpack import rfft, fftshift
import glob

import functools


import pickle
import random
import itertools


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
import gc

dask.set_options(get=dask.threaded.get)

########################hdf5 save and load #################################################################
def hd5save(df, name , overwrite ,verbose = False , pipelines= None ):
    #dataframe columns should all be arrays or bytestrings w ascii encoding

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

def applypipeline_to_series(series, pipeline, hyperparams):
    newseries = series.map( pipeline )
    if hyperparams['printResult']== True:
        print(newseries)

    newseries.index = series.index
    return newseries

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

def dataGen( fastas , fulldata = False):
    for fasta in fastas:
        fastaIter = SeqIO.parse(fasta, "fasta")
        for seq in fastaIter:
            if len(seq.seq)>0:
                if fulldata == False:
                    yield seq.seq
                else:
                    yield seq

###########################################################dataframe / dataset building ##########################################
def fastasToDF(fastas , DDF = None, verbose=False, ecodDB = False):
    regex = re.compile('[^a-zA-Z0-9]')
    regexAA = re.compile('[^ARDNCEQGHILKMFPSTWYV]')
    DFdict={}
    for fasta in fastas:
        if verbose == True:
            print(fasta)
        fastaIter = SeqIO.parse(fasta, "fasta")
        for seq in fastaIter:
            seqstr = regexAA.sub('', str(seq.seq))
            desc =str(seq.description)
            fastastr = '>'+desc+'\n'+seqstr+'\n'
            DFdict[desc] = { 'desc': desc.encode(), 'seq':seqstr, 'fasta': fastastr}
            if ecodDB == True:
                labels = ['ECOD uid','ECOD domain' , 'EOCD hierearchy string', 'ECOD pdb_range']
                for i,ecodb in enumerate(seq.description.split('|')[1:]):
                    DFdict[desc][labels[i]] = ecodb
    df = pd.DataFrame.from_dict(DFdict, orient = 'index')
    if verbose == True:
        print(df)
    
    DDF = dd.from_pandas(df , npartitions = 4*mp.cpu_count()  )
    
    if ecodDB == True:
        DDF.set_index('ECOD uid')
    else:
        DDF.set_index('desc')
    return DDF

def iter_sample_fast(iterator, samplesize):
    results = []
    # Fill in the first samplesize elements:
    for _ in range(samplesize):
        results.append(next(iterator))
    random.shuffle(results)  # Randomize their positions
    for i, v in enumerate(iterator, samplesize):
        r = random.randint(0, i)
        if r < samplesize:
            results[r] = v  # at a decreasing rate, replace random items

    if len(results) < samplesize:
        raise ValueError("Sample larger than population.")
    return results
########################################################################vis

from matplotlib import pyplot as plt

def plot_confusion_matrix(cm, classes, normalize=False, title='Confusion matrix', cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")
    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.show()



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





def serializeHash(minhash, hyperparams):
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

