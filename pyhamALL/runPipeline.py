
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

with open( config.working_dir + "speciestree.nwk" , 'r') as treefile:
    species_tree = treefile.read()
t = Tree(species_tree, format =1 )
print(t)
#rewrite IDs and reload tree
#
species_tree = pyham.utils.get_newick_string(working_dir + "speciestree_hack.nwk", type="nwk")

#open up OMA
h5file = open_file(omadir + 'OmaServer.h5', mode="r") 
#setup db objects
dbObj = db.Database(h5file)

omaIdObj = db.OmaIdMapper(dbObj)

#load Fams

#load orthoxml

#map the IDhack function to loaded orthoxml

#generate pyham objects

#generate minhash

#generate matrix rows

#comile LSH from serialized minhashes

#save lsh objects

#DONE


