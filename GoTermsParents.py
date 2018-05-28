from goatools import obo_parser
import h5py
import numpy as np
import tables

from utils import config_utils


def goterm2id(go_term_to_modif):
    id = int(go_term_to_modif.split(':')[1])
    return id


def _get_entry_gene_ontology(oma, entry):
    return oma.root.Annotations.GeneOntology.read_where('EntryNr == {}'.format(entry))


def _get_hog_members(hog_id, oma):
    """
    Gets all gene members from the hog
    :param hog_id: hog id
    :return: list of genes from the hog
    """
    iterator = _iter_hog_member(hog_id, oma)
    population = frozenset([x['EntryNr'] for x in iterator])
    return population


def _hog_lex_range(hog):
    """
    Decodes hog format
    :param hog: (bytes or string): hog id
    :return: hog_str: encoded hog id
    """
    hog_str = hog.decode() if isinstance(hog, bytes) else hog
    return hog_str.encode('ascii'), (hog_str[0:-1] + chr(1 + ord(hog_str[-1]))).encode('ascii')


def _iter_hog_member(hog_id, oma):
    """
    iterator over hog members / get genes
    :param hog_id: hog id
    :return: yields members of hog
    """
    hog_range = _hog_lex_range(hog_id)
    it = oma.root.Protein.Entries.where('({!r} <= OmaHOG) & (OmaHOG < {!r})'.format(*hog_range))
    for row in it:
        yield row.fetch_all_fields()


def iter_hog(oma):
    for hog in oma.root.HogLevel:
        yield hog[0]


def fam2hogid(fam_id):
    """
    Get hog id given fam
    :param fam_id: fam
    :return: hog id
    """
    hog_id = "HOG:" + (7-len(str(fam_id))) * '0' + str(fam_id)

    return hog_id


if __name__ == '__main__':

    obo_iterator = obo_parser.OBOReader(obo_file=config_utils.datadirLaurent + 'project/data/go.obo')
    obo_reader = obo_parser.GODag(obo_file=config_utils.datadirLaurent + 'project/data/go.obo')
    dt = h5py.special_dtype(vlen=np.dtype('int32'))
    omah5 = tables.open_file(config_utils.omadir + 'OmaServer.h5', mode='r')

    with h5py.File(config_utils.datadirLaurent + 'project/data/parents.h5', 'w', libver='latest') as h5_go_terms:
        # h5_go_terms.create_dataset('go_terms', (10000000,), dtype=dt)
        # dataset_go_terms_parents = h5_go_terms['go_terms']
        #
        # count = 0
        # for go_term in obo_iterator.__iter__():
        #     if go_term.namespace == 'biological_process':
        #
        #         try:
        #             go_term_read = obo_reader[go_term.id]
        #             go_term_parents = go_term_read.get_all_parents()
        #             go_term_parents_int = [goterm2id(go_term_read.id)] + [goterm2id(parent) for parent in go_term_parents]
        #             dataset_go_terms_parents[goterm2id(go_term_read.id)] = go_term_parents_int
        #
        #             count += 1
        #
        #             if count % 1000 == 0:
        #                 print('saving')
        #                 h5_go_terms.flush()
        #
        #         except KeyError:
        #             print('bug {}'.format(go_term.id))

        h5_go_terms.create_dataset('hog2genes', (1000000,), dtype=dt)
        dataset_hog2genes = h5_go_terms['hog2genes']

        h5_go_terms.create_dataset('genes2goterms', (biggest_gene,), dtype=dt)
        dataset_genes2goterms = h5_go_terms['genes2goterms']

        count = 0

        biggest_gene = 0

        for fam in iter_hog(omah5):
            genes2add = list(_get_hog_members(fam2hogid(fam), omah5))
            dataset_hog2genes[fam] = genes2add

            if biggest_gene < max(genes2add):
                biggest_gene = max(genes2add)
                print(biggest_gene)


            for gene in genes2add:
                dataset_genes2goterms[gene] = list(_get_entry_gene_ontology(omah5, gene))




            count =+ 1

            if count % 1000 == 0:
                print('saving')
                h5_go_terms.flush()

        print('the biggest gene is {}'.format(biggest_gene))

