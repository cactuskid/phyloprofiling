from goatools import obo_parser
import h5py
import numpy as np

from utils import config_utils


def goterm2id(go_term_to_modif):
    id = int(go_term_to_modif.split(':')[1])
    return id


if __name__ == '__main__':

    obo_iterator = obo_parser.OBOReader(obo_file=config_utils.datadirLaurent + 'project/data/go.obo')
    obo_reader = obo_parser.GODag(obo_file=config_utils.datadirLaurent + 'project/data/go.obo')
    dt = h5py.special_dtype(vlen=np.dtype('int32'))

    with h5py.File(config_utils.datadirLaurent + 'project/data/parents.h5', 'w', libver='latest') as h5_go_terms:
        h5_go_terms.create_dataset('go_terms', (10000000,), dtype=dt)
        dataset = h5_go_terms['go_terms']

        count = 0
        for go_term in obo_iterator.__iter__():
            if go_term.namespace == 'biological_process':

                try:
                    go_term_read = obo_reader[go_term.id]
                    go_term_parents = go_term_read.get_all_parents()
                    go_term_parents_int = [goterm2id(go_term_read.id)] + [goterm2id(parent) for parent in go_term_parents]
                    dataset[goterm2id(go_term_read.id)] = go_term_parents_int

                    count += 1

                    if count % 1000 == 0:
                        print('saving')
                        h5_go_terms.flush()

                except KeyError:
                    print('bug {}'.format(go_term.id))
