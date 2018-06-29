import profileGen
import h5py
import config
from datasketch import MinHashLSH , MinHashLSHForest
import datetime as dt

"""
generate a minhash lsh object based on a specific set of of evolutionary events


"""


datestr = "{:%B_%d_%Y_%H_%M}".format(dt.now())


lsh = MinHashLSH(threshold=0.7, num_perm=128)
forest = MinHashLSHForest(num_perm=128)
events = ['duplication', 'gain', 'loss', 'presence']

with open(config.datadir + datestr + 'lsh'.join(events)+'.pkl')
	with open(config.datadir + datestr + 'forest'.join(events)+'.pkl')
		for fam in h5hashes:
			try:
				minhash = profileGen.get_HOGhash(fam , h5hashes, events )
				lsh.insert(fam , minhash)
				forest.add(fam, minhash)
			except:
				print('error' + str(fam))
