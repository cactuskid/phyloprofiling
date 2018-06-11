import ujson as json
import h5py
import multiprocessing as mp
import pandas as pd
import random

from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import TermCounts

from utils import config_utils, hashutils
from validation import validation_semantic_similarity

from time import time

# data prep
go_terms_hdf5 = h5py.File(config_utils.datadirLaurent + 'project/data/parents.h5', 'r')
hogs2goterms = go_terms_hdf5['hog2goterms']


go = obo_parser.GODag(config_utils.datadirLaurent + 'project/data/go.obo')
associations = read_gaf(config_utils.datadirLaurent + 'project/data/gene_association.tair')

term_counts = TermCounts(go, associations)
goTermAnalysis = validation_semantic_similarity.Validation_semantic_similarity(go, term_counts, go_terms_hdf5)

def get_hogs_with_annotations(hogs):
    print('getting hogs')
    hogs_with_annotations = []
    hogs_without_annotations = []
    for fam, goterms in enumerate(hogs):

        try:
            obj = json.loads(goterms)
            if obj and type(obj) is dict and len(obj) > 0:
                hogs_with_annotations.append(fam)
            else:
                # print(obj)
                hogs_without_annotations.append(fam)
        except ValueError:
            hogs_without_annotations.append(fam)

        if len(hogs_without_annotations) > 100000:
            break
        # if goterms:
        #     if json.loads(goterms):
        #         hogs_with_annotations.append(hashutils.fam2hogid(fam))

    print('hogs without annotations {}'.format(len(hogs_without_annotations)))
    print('hogs with annotations {}'.format(len(hogs_with_annotations)))
    return hogs_with_annotations


start_time = time()
hogs_w_annotations = get_hogs_with_annotations(hogs2goterms)
print('time to get hogs: {}'.format(time()-start_time))

number_hogs_with_annotations = len(hogs_w_annotations)
small_hogs_with_annotations = random.sample(hogs_w_annotations, 1000)
print(small_hogs_with_annotations)


def make_dict(fam):
    try:
        print(fam)
        returnDict = {fam: goTermAnalysis.semantic_similarity_score(hashutils.fam2hogid(fam), hashutils.fam2hogid(fam))}
    except:
        print(" *** {} ***".format(fam))
        returnDict = {fam: -1}
    return returnDict

start_time = time()
if __name__ == '__main__':

    pool = mp.Pool()
    results = pool.map_async(make_dict, small_hogs_with_annotations, 10*mp.cpu_count())

    results = results.get()

    global_dict = {}
    for r in results:
        global_dict.update(r)

    dt = pd.DataFrame.from_dict(global_dict, orient='index')
    dt.to_csv('hogvshog.csv', sep='\t')

print(time()-start_time)
