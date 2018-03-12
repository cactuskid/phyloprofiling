

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

#species_tree = pyham.utils.get_newick_string(working_dir + "speciestree.nwk", type="nwk")
with open( config.working_dir + "speciestree.nwk" , 'r') as treefile:
    species_tree = treefile.read()


h5file = open_file(omadir + 'OmaServer.h5', mode="r") 
#setup db objects
dbObj = db.Database(h5file)
omaIdObj = db.OmaIdMapper(dbObj)


t = Tree(species_tree, format =1 )
print(t)
species_tree = pyham.utils.get_newick_string(working_dir + "speciestree_hack.nwk", type="nwk")



