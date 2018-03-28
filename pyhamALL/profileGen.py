#use dask and pyham to create a big sparse matrix for each type of evolutionary event
import pyham
import dask
import datasketch
import ete3
import sparse
import itertools
from datasketch import LeanMinHash, MinHash, MinHashLSH

def generateTaxaIndex(newick):
	'''
	Generates taxa index
	Generates an index for the global taxonomic tree for all OMA
	Args:
		newick: species tree in newick format
	Returns:
		taxaIndex: dictionary key: node name (species name); value: index
		taxaIndexReverse: dictonary key: index: value: species name
	'''
	t =ete3.Tree(newick)
	taxaIndex = {}
	for i,node in enumerate(t.traverse):
		taxaIndexReverse[i] = node.name
		taxaIndex[node.name] = i
	return taxaIndex, taxaIndexReverse



def Tree2Hashes(eteobj):
	#turn each tree into a minhash object
	#serialize and store as array
	eventdict = { 'presence':[] , 'gain':[] , 'loss':[] , 'duplication':[]}	



	eventdict = { 'presence':[] , 'gain':[] , 'loss':[] , 'duplication':[]}
	

	for node in eteobj.traverse():
	# traverse() returns an iterator to traverse the tree structure
	# strategy:"levelorder" by default; nodes are visited in order from root to leaves
	# it return treeNode instances
		node.

	hashes= []

	for array in eventdict:
		#generate minhash
		# why set ? 
		eventdict[array] = set(eventdict[array])
		m1 = MinHash(num_perm=128)
		for element in eventdict[array]:
			m1.update(element)
		m1 = LeanMinHash(m1)
		buf = bytearray(m1.bytesize())
		m1.serialize(buf)
		hashes.apped([m1.serialize()])
	hashmat = np.vstack(hashes)
	return hashmat


def Tree2mat(eteobj, taxaIndex):
	#use partials to configure the taxa index
	#turn each tree into a sparse matrix with 4 rows
	rowdict={ 'presence':0 , 'gain':1 , 'loss':2 , 'duplication':3}
	matrix = sparse((len(rowdict), len(taxaIndex) ))
	
	for node in tree.traverse():
		
		matrix[ rowdict[] , taxaIndex[] ] = 1


def MatToLSH(index , hashmat , LSH , rownum = None):
	#lsh forest single core?
	if rownum == None:
		for row in range(hashmat.shape[0]):
			#deserialize byte array and insert it into lsh
			m1 = LeanMinHash.deserialize(hashmat[row,:])
			LSH.insert(index[row], m1)
	else:
		#build of an lsh with only a few of the hash rows
		for row in rownum:
			m1 = LeanMinHash.deserialize(hashmat[row,:])
			LSH.insert(index[row], m1)
	#use itertools to build all 15 combos of hash signatures using serialized leanminhashes

	



