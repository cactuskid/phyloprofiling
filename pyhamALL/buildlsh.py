import profileGen
import h5py
import config
from datasketch import MinHashLSH , MinHashLSHForest

lsh = MinHashLSH(threshold=0.7, num_perm=128)
forest = MinHashLSHForest(num_perm=128)
for fam in h5hashes:
	try:
		minhash = get_HOGhash(fam , h5hashes, events = ['duplication', 'gain', 'loss', 'presence'])
		lsh.insert(fam , minhash)
		forest.add(fam, minhash)
	except:
		print('error' + str(fam))
		