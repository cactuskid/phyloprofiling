
# coding: utf-8

# The purpose of this notebook is to address some of the issues with running pyham in its current version. It includes work-arounds and tips to run pyham with the current version of the OMA database. At the end shows a simple analysis on how to use pyham to find the number of duplicated or lost genes for a given HOG.
# 
# by Natasha Glover, last updated 11 Jan 2018

# # Get setup

# In[2]:


import pyham
from pyoma.browser import db 
import numpy as np
from tables import *
import re
from ete3 import Tree
import pickle
import mpwtimeout as multi
# In[3]:
import tempfile
import functools



working_dir = "./"
datadir = '/scratch/cluster/monthly/dmoi/dmoiProfiling/'
omadir = '/scratch/ul/projects/cdessimo/oma-browser/All.Dec2017/data/'

buildtestdataset = True


#species_tree = pyham.utils.get_newick_string(working_dir + "speciestree.nwk", type="nwk")
with open( working_dir + "speciestree.nwk" , 'r') as treefile:
    species_tree = treefile.read()

#species_tree = ete3. (datadir + "speciestree.nwk", type="nwk")

# Here are some useful functions:

# In[5]:


def replace_characters(string):
    string = string.replace(".", "_")
    string = string.replace(" ", "_")
    string = string.replace("(", "")
    string = string.replace(")", "")
    string = string.replace(":", "-")
    return(string)

def fix_species_tree(species_tree):
    '''replaces characters which mess up the newick species tree'''
    replacement_dic = {}
    species_and_taxa = re.findall(r'"([^"]*)"', species_tree) 
    
    for old_name in species_and_taxa:
        new_name = replace_characters(old_name)
        replacement_dic[old_name] = new_name

    #sort by length of value so that long items get replaced first
    old_names_list = list(replacement_dic.keys())
    sorted_old_names_list = sorted(old_names_list , key = len, reverse=True)

    for name in sorted_old_names_list:
        species_tree = species_tree.replace( name, replacement_dic[name])

    species_tree = species_tree.replace("\"", "")
    species_tree = species_tree.replace("\n", "")
    
    return(species_tree)

def get_species_sciname(uniprotspeciescode):
    sciname = genome_df[genome_df['UniProtSpeciesCode']==uniprotspeciescode.encode("utf-8")]['SciName'].item().decode("utf-8") 
    return(sciname)


# In[6]:

#replace characters species tree

species_tree = fix_species_tree(species_tree)
with open( working_dir + 'speciestree_hack.nwk' , 'w') as outTree:
    outTree.write(species_tree)


t = Tree(species_tree, format =1 )
print(t)
species_tree = pyham.utils.get_newick_string(working_dir + "speciestree_hack.nwk", type="nwk")


h5file = open_file(omadir + 'OmaServer.h5', mode="r") 

#setup db objects
dbObj = db.Database(h5file)
omaIdObj = db.OmaIdMapper(dbObj)

#get a list of all the species which are identified as their 5-letter uniprot species code in the species tree
uniprot_species = []
uniprot_species = re.findall(r'\b[A-Z]{5}\b', species_tree)
uniprot_species.append("STAA3")
uniprot_species.append("ECO57")
uniprot_species.append("BUCAI")
uniprot_species.append("CHLPN")

print("UniProt 5-letter species codes which have replaced their scientific names in the species tree: "+str(uniprot_species))

#make a dictionary with the scientific name and the uniprot species code
#this replacement dic is for later to replace key with value in orthoxml

replacement_dic ={}
for species in uniprot_species:
    try:
        replacement_dic[replace_characters(omaIdObj.genome_from_UniProtCode(species)[5].decode())] = species

    except:
        pass
    
print("\nreplacement_dic: "+ str(replacement_dic))


# # Getting the orthoxml file

# There are two ways that I know of to get an orthoxml file. 
# 
# 1) using pyoma to extract the orthoxml for a particular HOG from the HDF5 file, or
# 
# 2) downloading the whole orthoxmlfile for the HOGs from the omabrowser: https://omabrowser.org/oma/current/
# 
# Since option 2 will yield a really HUGE file which will take several hours to create the ham object, I recommend and will focus on option 1 for this tutorial.

# Not so fast, you're not done with this orthoxml file. Still have to make a few modifications-- replacing all the special characters and species names from the previous section.

# In[23]:


def convert_orthoxml_ids(myinfile, myoutfile, replacement_dic):
    '''Takes an orthoxml file as input, along with the replacement_dic, where the keys are scientific names (with 
    special characters already replaced) and values are the new name which matches the species tree.
    Replaces the special characters and old scientific names with the new ones.
    Returns the number of genes in the orthoxml file.'''
    
    mylist = []
    outfile = open(myoutfile , "w")
    count = 0 
    
    with open(myinfile) as infile:
        for line in infile:
            searchObj = re.search( r'.*<species name=\"(.*)\" NCBITaxId.*', line)
            searchObj2 = re.search(r'.*<property name=\"TaxRange\" value=\"(.*)\"\/>', line)
            searchObj3 = re.search(r'.*<gene id=.*', line)
            if searchObj:
                old_name = searchObj.group(1)
                new_name = replace_characters(old_name)
                line = line.replace(old_name, new_name)

                for key, value in replacement_dic.items():
                    if new_name == key:
                        line = line.replace(key, value)

            if searchObj2:
                old_name = searchObj2.group(1)
                new_name = replace_characters(old_name)
                line = line.replace(old_name, new_name)

                for key, value in replacement_dic.items():
                    if new_name == key:
                        line = line.replace(key, value)   
                        
            if searchObj3:
                count = count + 1

            outfile.write(line)      
        
    outfile.close()
    return(count)

#pyham takes a file rather than a string, so save it as a local file
def retham(fam, l, dbObj, species_tree, datadir , replacement_dic):
    dbObj.get_orthoxml(fam)
    #make ham analysis object
    with tempfile.NamedTemporaryFile(dir = datadir ) as temp:
        l.acquire()
        ortho = dbObj.get_orthoxml(fam)
        l.release()
        temp.write( ortho )

        with tempfile.NamedTemporaryFile(dir = datadir ) as temp2:
            try:
                index = 'HOG:'.join(['0']*(6-len(str(fam))) + fam )
                l.acquire()
                nb_genes = convert_orthoxml_ids(myinfile = datadir+temp.name , 
                         myoutfile =  datadir + temp2.name ,
                         replacement_dic = replacement_dic)
                l.release()
                hamObj = pyham.Ham( species_tree, datadir + temp2.name )
                print(str(index)+':ham done')
                return {index: hamObj}
            except:
                print ('pyham error')
                print (str(fam))



def retham_testdataset(fam,  dbObj, species_tree, testdir , replacement_dic):
    #make ham analysis object
    index = str(fam)
    with open( testdir + index +'.orthoxml' , 'w' ) as outfile:
        ortho = dbObj.get_orthoxml(fam)
        print(ortho)
        outfile.write( str(ortho) )
    
    nb_genes = convert_orthoxml_ids(myinfile = testdir + index +'.orthoxml'  , 
             myoutfile =  testdir + index +'_IDhack.orthoxml'  ,
             replacement_dic = replacement_dic)
    
    hamObj = pyham.Ham( species_tree, testdir + index +'_IDhack.orthoxml'  , use_internal_name= True , format = 1)
    print(str(index)+':ham done')
    return {index: hamObj}

if buildtestdataset == True:
    testdir = './test/'
    hamdict={}
    for row in h5file.root.OrthoXML.Index[0:100]:
        hamdict.update(retham_testdataset(row[0], dbObj , species_tree , testdir , replacement_dic ))
    with open( testdir + 'hamdict.pkl' , 'wb' ) as handle:
        pickle.dump(hamdict , handle , -1 )

else:
    def yeildfams():
        for row in h5file.root.OrthoXML.Index:
            yield row[0]
    retHamMP = functools.partial( retham , dbObj=dbObj , species_tree= species_tree,  datadir = datadir , replacement_dic = replacement_dic )
    multi.mp_with_timeout(nworkers= 10, nupdaters = 1, startobject ={} , saveobject= datadir + 'hams.hdf5'  , 
        datagenerator= yeildfams()  , workerfunction = retHamMP, updaterfunction=multi.daskupdater ,updateobjfunction =multi.updatefunction_dict , timeout = 60, saveinterval = 600  )



