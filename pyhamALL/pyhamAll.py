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

# with open( config.working_dir + "speciestree.nwk" , 'r') as treefile:
#     species_tree = treefile.read()
# t = Tree(species_tree, format =1 )
# print(t)
# #rewrite IDs and reload tree
# #
# species_tree = pyham.utils.get_newick_string(working_dir + "speciestree_hack.nwk", type="nwk")


#open up OMA
h5file = open_file(config.omadirLaurent + 'OmaServer.h5', mode="r") 
#setup db objects
dbObj = db.Database(h5file)
omaIdObj = db.OmaIdMapper(dbObj)

# corrects species tree and replacement dictionary for orthoXML files
species_tree, replacement_dic = formatOrtho.fix_species_tree("speciestree.nwk", omaIdObj)

#temporary, 
hashmat_list = [] 
mat_list = []

# init LSH

#load Fam
for fam in pyhamPipeline.yieldFamilies(h5file):
	try:
	# generate pyham object
	ham_fam = pyhamPipeline.get_ham(fam, dbObj, species_tree, replacement_dic)

	# generate tree profile
	treemap_fam = pyhamPipeline.pyhamtoTree(ham_fam, fam)

	# generate matrix of hash
	hashmat = profileGen.Tree2Hashes(treemap_fam)
	hashmat_list.append(hashmat)
	# generate taxa index
	taxaIndex, taxaIndexReverse = profileGen.generateTaxaIndex(species_tree)

	# generate matrix of 1 and 0 for each biological event
	mat = profileGen.Tree2mat(treemap_fam, taxaIndex)
	mat_list.append(mat)
	print(fam)
except:
	pass

#load orthoxml

#map the IDhack function to loaded orthoxml

#generate pyham objects

#generate minhash

#generate matrix rows

#compile LSH from serialized minhashes

#save lsh objects

#DONE


