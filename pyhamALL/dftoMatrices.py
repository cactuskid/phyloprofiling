#use dask and pyham to create a big sparse matrix for each type of evolutionary event
import pyham
import dask
import datasketch
import ete3


def generateTaxaIndex(newick):
	#use this to generate an index for the global taxonomic tree for all of OMA
	t =ete3.Tree( newick)
	taxaIndex = {}
	for i,n in enumerate(t.traverse):
		taxaIndexReverse[i] = node.name
		taxaIndex[node.name] = i
	return taxaIndex, taxaIndexReverse

def pyhamtoTree(hamOBJ):
	#use pyham to get all evolutionary events from a pyham object
	#turn into an array and a hash
	tp = pyham.TreeProfile(hamOBJ)
	tree = tp.compute_tree_profile_full()
	return tree

def Tree2Hashes(eteobj):
	#turn each tree into a minhash object
	presence = []
	gain = []
	loss = []
	duplication = []

	for node in tree.traverse():
		node.


	return presence , gain, loss, duplication

def Tree2mat(eteobj, taxaIndex):
	#turn each tree into a matrix row
	for node in tree.traverse():
		node.




