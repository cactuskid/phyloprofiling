import ujson as json
import h5py
import multiprocessing as mp
import pandas as pd

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
    for fam, goterms in enumerate(hogs):
        if goterms:
            if json.loads(goterms):
                hogs_with_annotations.append(hashutils.fam2hogid(fam))


    print(len(hogs_with_annotations))
    return hogs_with_annotations
start_time = time()
hogs_w_annotations = get_hogs_with_annotations(hogs2goterms)
print('time to get hogs: {}'.format(time()-start_time))

def fill_dict(dict, hog):
    dict[hog] = goTermAnalysis.semantic_similarity_score(hog, hog)


#for hog in hogs_w_annotations:
#    hogs_semantic_score[hog] = goTermAnalysis.semantic_similarity_score(hashutils.fam2hogid(hog), hashutils.fam2hogid(hog))

start_time = time()
if __name__ == '__main__':

    # mgr = mp.Manager()
    # hogs_semantic_score = mgr.dict()
    # job = [mp.Process(target=fill_dict, args=(hogs_semantic_score, hog)) for hog in hogs_w_annotations]
    # _ = [p.start() for p in job]
    # _ = [p.join() for p in job]
    hogs_semantic_score = {}

    for hog in hogs_w_annotations:
        fill_dict(hogs_semantic_score,hog)


    print(hogs_semantic_score)

    dt = pd.DataFrame.from_dict(hogs_semantic_score, orient='index')
    dt.to_csv('hogvshog.csv', sep='\t')

print(time()-start_time)