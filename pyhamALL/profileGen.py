#use dask and pyham to create a big sparse matrix for each type of evolutionary event
import pyham
import dask
import ete3
import sparse
from scipy.sparse import csr_matrix,find 
import itertools
import datasketch
import numpy as np
import itertools

def generateTaxaIndex(species_tree):
	'''
	Generates an index for the global taxonomic tree for all OMA
	Args:
		species_tree: species tree in newick format
	Returns:
		taxaIndex: dictionary key: node name (species name); value: index
		taxaIndexReverse: dictonary key: index: value: species name
	'''
	t =ete3.Tree(species_tree, format=1)
	taxaIndex = {}
	taxaIndexReverse = {}
	for i,node in enumerate(t.traverse()):
		taxaIndexReverse[i] = node.name
		taxaIndex[node.name] = i
	return taxaIndex, taxaIndexReverse

def FamList2RowsOnTheFly(listfam):
	pass

def FamList2RowsH5(h5file, listfam):
	pass



def jaccard_cutoff(fams, scores, cutoff):
    return fams[ np.where(scores > cutoff)]


def get_HOGhash(fam , h5hashes. events = ['duplication', 'gain', 'loss', 'presence']):
	#get hash of desired events
    for event in events:
            buf = h5hashes[event][fam,:]
            if queryminhash is None:
                query_minhash = LeanMinHash.deserialize(buf)                
                minhash1 = MinHash(seed=quers_minhash.seed, hashvalues=query_minhash.hashvalues)
            else:
                query_minhash =  LeanMinHash.deserialize(buf)
                minhash2 = MinHash(seed=quers_minhash.seed, hashvalues=query_minhash.hashvalues)
                minhash1.merge(minhash2)
    return minhash1

def get_hashdict(fams, h5hashes, events = ['duplication', 'gain', 'loss', 'presence']):
	hashdict = {}
	for fam in fams:
		hashdict[fam] = get_HOGhash( fam, h5hashes, events)
	return hashdict
	
def jaccard_rank(query,hashdict):
    hashlist = np.asarray(list(hashdict.values()))
    fams = np.asarray(list(hashdict.keys()))
    scores = [ queryhash.jaccrad(resulthash) for resulthash in hashlist ]
    index = np.argsort(scores)
    for array in [hashlist, fams , index]:
        array = array[index]
    return fams,scores


def Tree2Hashes(treemap, fam=None, LSH=None , l = None):
	#turn each tree into a minhash object
	#serialize and store as array
	eventdict = { 'presence':[] , 'gain':[] , 'loss':[] , 'duplication':[]}	
	

	for node in treemap.traverse():
	# traverse() returns an iterator to traverse the tree structure
	# strategy:"levelorder" by default; nodes are visited in order from root to leaves
	# it return treeNode instances
		if not node.is_root():
			if node.nbr_genes >0:
				eventdict['presence'].append('P'+node.name)
			if node.dupl > 0:
				eventdict['duplication'].append('D'+node.name)
			if node.lost > 0:
				eventdict['loss'].append('L'+node.name)
		else:
			eventdict['gain'].append('G'+node.name)

	hashes = {}

	hashesDict = {}

	lminHashDict = {}

	for array in eventdict:
		eventdict[array] = set(eventdict[array])

		minHash = datasketch.MinHash(num_perm=128)

		for element in eventdict[array]:
			
			minHash.update(element.encode())
			

		hashesDict[array] = minHash

		lminHash = datasketch.LeanMinHash(minHash)

		lminHashName = str(fam)+'-'+array

		lminHashDict[lminHashName] = lminHash
		if LSH is not None: 
			LSH.insert(lminHashName, lminHash)

		
		buf = bytearray(lminHash.bytesize())
		lminHash.serialize(buf)

		hashes[array] = buf


	#hashmat = np.hstack(buffers)

	for j in range(1,len(eventdict.keys())):
		for i in itertools.combinations(eventdict.keys(), j+1):
			combName = str(fam)
			minHash = datasketch.MinHash(num_perm=128)
			for array in i:
				combName += '-'+array
				minHash.merge(hashesDict[array])

			lminHashDict[combName] = minHash
			lminHash = datasketch.LeanMinHash(minHash)	
		if LSH:
			LSH.insert( combName , lminHash)

	return hashes , lminHashDict

def DFTree2Hashes(row):
	#turn each tree into a minhash object
	#serialize and store as array
	
	fam, treemap = row.tolist()

	eventdict = { 'presence':[] , 'gain':[] , 'loss':[] , 'duplication':[]}	
	if treemap is not None:
		for node in treemap.traverse():
		# traverse() returns an iterator to traverse the tree structure
		# strategy:"levelorder" by default; nodes are visited in order from root to leaves
		# it return treeNode instances
			if not node.is_root():
				if node.nbr_genes >0:
					eventdict['presence'].append('P'+node.name)
				if node.dupl > 0:
					eventdict['duplication'].append('D'+node.name)
				if node.lost > 0:
					eventdict['loss'].append('L'+node.name)
			else:
				eventdict['gain'].append('G'+node.name)

		hashes = {}

		hashesDict = {}

		lminHashDict = {}

		for array in eventdict:
			eventdict[array] = set(eventdict[array])

			minHash = datasketch.MinHash(num_perm=128)

			for element in eventdict[array]:
				
				minHash.update(element.encode())
				

			hashesDict[array] = minHash

			lminHash = datasketch.LeanMinHash(minHash)

			lminHashName = str(fam)+'-'+array

			lminHashDict[lminHashName] = lminHash
			
			
			buf = bytearray(lminHash.bytesize())
			lminHash.serialize(buf)

			hashes[array] = buf


		#hashmat = np.hstack(buffers)

		for j in range(1,len(eventdict.keys())):
			for i in itertools.combinations(eventdict.keys(), j+1):
				combName = str(fam)
				minHash = datasketch.MinHash(num_perm=128)
				for array in i:
					combName += '-'+array
					minHash.merge(hashesDict[array])

				lminHashDict[combName] = minHash
				lminHash = datasketch.LeanMinHash(minHash)	
			
		return {'hashes':hashes , 'dict': lminHashDict}
	else:
		return None

def Tree2mat( treemap, taxaIndex, verbose = False):
	'''
	Turn each tree into a sparse matrix with 4 rows.

	Args:
		treemap : tree profile object from pyHam; contains biological events
		taxaIndex : index of gloval taxonomic tree from OMA

	Returns :
		profile_matrix : matrix of size numberOfBiologicalEvents times taxaIndex containing when the given biological event is present in the given species
	'''
	#use partials to configure the taxa index
	if treemap is not None:
		hogmat = csr_matrix( (1 , 4*len(taxaIndex)  ) )
		columndict={ 'presence':0
		 , 'gain': 1 
		 , 'loss': 2
		 , 'duplication':3
		}

		
		for node in treemap.traverse():
		# traverse() returns an iterator to traverse the tree structure
		# strategy:"levelorder" by default; nodes are visited in order from root to leaves
		# it return treeNode instances
			if not node.is_root():
				# for presence, loss, and duplication, set 1 if > 0           
				if node.nbr_genes > 0:
					pad = columndict['presence'] *len(taxaIndex)
					hogmat[ 0 , pad + taxaIndex[node.name]] = node.nbr_genes 
				if node.lost > 0:
					pad = columndict['loss']*len(taxaIndex)
				
					hogmat[0 , pad + taxaIndex[node.name]] = 1
				if node.dupl > 0:
					pad = columndict['duplication']*len(taxaIndex)
					hogmat[0 , pad + taxaIndex[node.name]] = 1
			else:
				# gain is only for root; impossible to "gain" a gene several times
				pad = columndict['gain']*len(taxaIndex)
				hogmat[0 , pad + taxaIndex[node.name]] = 1
		
		if verbose == True:
			print(hogmat)

		return hogmat

	else:
		return None


