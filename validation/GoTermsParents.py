from goatools import obo_parser
import h5py
import numpy as np

from utils import config_utils


obo_reader = obo_parser.OBOReader(obo_file=config_utils.datadirLaurent + 'project/data/go.obo')
dt = h5py.special_dtype(vlen=np.dtype('int32'))

with h5py.File('saving path', 'w', libver='latest') as h5_go_terms:
    h5_go_terms.create_dataset('go_terms', (10000000,), dtype=dt)
    dataset = h5_go_terms['go_terms']

    i for enumerate go terms
    dataset[54] = find parents for go:000054

    if i % 1000 == 0:
        every 1000h5_go_terms.flush()


# loop over go terms in obo file (take only biological process)

# for each go term -> get all the parents