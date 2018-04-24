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



h5file_save = h5py.File(config.datadir + 'hogProfiles', 'a')

dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']

# TODO change chucksize for full run
chunksize = 3

for dataset_name in dataset_names:
	if dataset_name not in list(h5file_save.keys()):
		dataset = h5file_save.create_dataset(dataset_name, (0,0), maxshape=(None, None), dtype = 'int32')

dsets = {}
for dataset_name in list(h5file_save.keys()):
	dsets[dataset_name] = h5file_save[dataset_name]


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

import pyhamPipeline
import profileGen
import format_files

from datasketch import MinHashLSH

import h5py
import h5sparse





parallel = False
#open up OMA
h5file = open_file(config.omadir + 'OmaServer.h5', mode="r") 
#setup db objects
dbObj = db.Database(h5file)
omaIdObj = db.OmaIdMapper(dbObj)
# corrects species tree and replacement dictionary for orthoXML files
species_tree, replacement_dic = formatOrtho.fix_species_tree(h5file, omaIdObj)
taxa_index = profileGen.generateTaxaIndex(species_tree)

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
	hashmat_list = [] 
	mat_list = []


			for dset in dsets:
				if dset == 'fam':
					dsets[dset][i,:] = fam
				else:
					dsets[dset][i,:] = hashesDic[dset]


h5file = h5py.File(config.omadirLaurent + 'data1', 'a')
h5fmatrix =  h5sparse.File(config.datadir + "matrix.h5")


dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']
list(h5file.keys())
for dataset_name in dataset_names:
    if dataset_name not in list(h5file.keys()):
        dataset = h5file.create_dataset(dataset_name, (0,0), maxshape=(None, None), dtype = 'int32')

#open up OMA
h5file = open_file(config.omadirLaurent + 'OmaServer.h5', mode="r") 
#setup db objects
dbObj = db.Database(h5file)
omaIdObj = db.OmaIdMapper(dbObj)

# corrects species tree and replacement dictionary for orthoXML files
dic, tree = format_files.create_species_tree(h5file, omaIdObj)
#set up tree mapping dictionary
taxaIndex = profileGen.generateTaxaIndex(tree)

#initialize matrix h5 file
for dataset_name in dataset_names:
    treemap_fam = pyhamPipeline.get_hamTree(0, dbObj, tree, dic)
    matrixRow = profileGen.Tree2mat(treemap_fam, taxaIndex)
    if dataset_name != 'fam':
        dataset = h5fmatrix.create_dataset(dataset_name, data=matrixRow)
    else:
        dataset = h5fmatrix.create_dataset(dataset_name, (0,0), maxshape=(None, None), dtype = 'int32')
        

dsets = {}
matrixdsets = {}
for dataset_name in list(h5file.keys()):
    dsets[dataset_name] = h5file[dataset_name]
    matrixdsets[dataset_name] = h5fmatrix[dataset_name]




chunksize = 3
#temporary, 
mat_list = []
dico_list = []
lsh = MinHashLSH()

#load Fam
for i,fam in enumerate(pyhamPipeline.yieldFamilies(h5file)):
    
    if fam > 0 and fam < 5:
        print(fam)
    
    #generate treemap profile
    treemap_fam = pyhamPipeline.get_hamTree(fam, dbObj, tree, dic)
    # generate matrix of hash
    hashesDic = profileGen.Tree2Hashes(treemap_fam, fam, lsh)
    matrixRow = profileGen.Tree2mat(treemap_fam, taxaIndex)


    if i == 0:
        for dataset in dsets:
            if dataset == 'fam':
                dsets[dataset].resize((chunksize, 1 ))
            else: 
                dsets[dataset].resize((chunksize, np.asarray(hashesDic[dataset]).shape[0] ))
                



    if i % chunksize == 0 and i != 0:
        for dataset in dsets:
            if dataset == 'fam':
                dsets[dataset].resize((i+chunksize, 1 ))
            else: 
                dsets[dataset].resize((i+chunksize, np.asarray(hashesDic[dataset]).shape[0] ))

    for dset in dsets:
        if dset == 'fam':
            dsets[dset][i,:] = fam
        else:
            dsets[dset][i,:] = hashesDic[dset]


	h5file.close()
#         dico_list.append(dico)
#         # generate taxa index
#         taxaIndex, taxaIndexReverse = profileGen.generateTaxaIndex(tree)

#         # generate matrix of 1 and 0 for each biological event
#         mat = profileGen.Tree2mat(treemap_fam, taxaIndex)
#         mat_list.append(mat)

#     except:
#         print("breaking at {}".format(fam))
        

#load orthoxml

#map the IDhack function to loaded orthoxml

#generate pyham objects

#generate minhash

#generate matrix rows

#compile LSH from serialized minhashes

#save lsh objects

#DONE


