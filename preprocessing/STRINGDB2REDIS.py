import redis
import gc
from utils import config_utils


UNIMAP = True
STRINGMAP = True


if STRINGMAP:
    r1 = redis.StrictRedis(host='localhost', port=6379, db=1)
    # sort the IDs alphanumerically.
    # protein1 protein2 neighborhood fusion cooccurence coexpression
    # experimental database textmining combined_score
    # 394.NGR_c00010 394.NGR_c33930 0 0 165 0 0 0 145 255
    # 37200000
    # save file line...
    refs = ['neighborhood', 'fusion', 'cooccurence', 'coexpression', 'experimental', 'database', 'textmining', 'combined_score']
    start_line = 0
    with open(config.string_interactors +'/stringdata/protein.links.detailed.v10.5.txt', 'r') as stringAll:
        for i, line in enumerate(stringAll):
            if i > start_line:
                words = line.split()
                IDS = ''.join(sorted([words[0], words[1]]))
                r1.set(IDS, stringAll.tell())
            if i % 1000000 == 0:
                print(i)
				if i % 10000000 == 0:
					gc.collect()
