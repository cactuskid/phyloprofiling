#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 17:13:21 2018

@author: laurent
"""
import phyloValidationGoTerm

import tables
import numpy as np

from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import TermCounts 

from pyoma.browser import db

# load all the files (oma database, go dag file, association gaf file)
# set oma database
# contains genes and related GO terms
omaDataLocation = "/home/laurent/Documents/phylo/phyloprofiling_working/1stdraft/data/OmaServer.h5"
h5file = tables.open_file(omaDataLocation, mode='r')

dbObj = db.Database(h5file)
omaIdObj = db.OmaIdMapper(dbObj)
#load gene ontology DAG
#use guaranteed acyclic basic GO
go = obo_parser.GODag('/home/laurent/Documents/phylo/phyloprofiling_working/1stdraft/data/go.obo')
associations = read_gaf('/home/laurent/Documents/phylo/phyloprofiling_working/1stdraft/data/gene_association.tair.gz')

# Get the counts of each GO term.
termcounts = TermCounts(go, associations)

# random hogs generation for testing
# test the analysis, take hogs form drosophila melanogaster
drosoHogList = []
for hog in dbObj.get_hdf5_handle().get_node('/HogLevel').where('Level == b"Drosophila melanogaster"'):
    drosoHogList.append(hog[1])

# --------------------------------------------------------
# ---------------------- With Genes ----------------------
# --------------------------------------------------------
  
# create object for validation; set files
validation_hogs = phyloValidationGoTerm.SemanticSimilarityAnalysis(go, h5file, termcounts)

# take a few hogs with annotations
hogsList =[]

import time

t0 = time.time()
for hog in drosoHogList:
    if len(hogsList) > 50:
        break
    if len(validation_hogs.get_go_terms(hog)) > 15:
        hogsList.append(hog)
print(time.time()-t0)

validation_hogs.semantic_similarity_score(hogsList[3],hogsList[3])