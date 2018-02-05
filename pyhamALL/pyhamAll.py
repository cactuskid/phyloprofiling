
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

#here I use the Mar2017 release of the OMA database
h5file = open_file(omadir + 'OmaServer.h5', mode="r") 

#setup db objects
dbObj = db.Database(h5file)
omaIdObj = db.OmaIdMapper(dbObj)
#taxObj = db.Taxonomy(np.array(h5file.root.Taxonomy[:]))   


species_tree = pyham.utils.get_newick_string(working_dir + "speciestree.nwk", type="nwk")


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

print(species_tree)


# # Figuring out which species names have been replaced

# You can compare the old version and the new version of the species tree to see the special characters that have been replaced. **HOWEVER, ** there is another problem. Certain scientific names in the species tree have been replaced by the 5-letter species code. This is the case when an internal node has the same name as a species (leaf). An example is E.coli. If you use Ctrl+F to find "Escherichia_coli_strain_K12" in the above species tree, you will find 3 matches-- 2 of these are species, and the last one is the internal node id. Hence, the actual Escherichia_coli_strain_K12 species itself has been replaced by the 5-letter UniProt code, which is **ECOLI**. Since some of the species names have been replaced in the newick species tree, they must be replaced in the orthoxml. 
# 
# The following code gets all these UniProt species codes found in the newick tree, and collects them in a list (uniprot_species). Then, it creates a dictionary where the key is the string to be replaced in the orthoxml and the value is the uniprot species code.

# In[7]:


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
        try:

            nb_genes = convert_orthoxml_ids(myinfile = temp.name , 
                     myoutfile = datadir + str(fam)+'ALLhogs_IDhack.orthoxml',
                     replacement_dic = replacement_dic)
            hamObj = pyham.Ham( species_tree, datadir + str(fam)+'ALLhogs_IDhack.orthoxml' , use_internal_name=True)
            print(str(fam)+':ham done')
            
            return {fam: hamObj}
        except:
            print ('pyham error')
            print (str(fam))

print(h5file.root.OrthoXML.Index[0:10])

def yeildfams():
    for row in h5file.root.OrthoXML.Index:
        yield row[0]

retHamMP = functools.partial( retham , dbObj=dbObj , species_tree= species_tree,  datadir = datadir , replacement_dic = replacement_dic )
multi.mp_with_timeout(nworkers= 10, nupdaters = 1, startobject ={} , saveobject=datadir + 'hams.pkl'  , 
    datagenerator= yeildfams()  , workerfunction = retHamMP, updatefunction =multi.updatefunction , timeout = 60, saveinterval = 600  )



