from goatools import obo_parser
import h5py
import numpy as np
import tables
import ujson as json
from pyoma.browser import db
import redis
import gc
from time import time

from utils import config_utils
from utils import preprocess_config
import StringRedisTOOLS

def yield_hogs_with_annotations(annotation_dataset):
    for fam, annotations in enumerate(hogs):
        try:
            obj = json.loads(annotations)
            if obj and type(obj) is dict and len(obj) > 0:
                yield {fam : obj}
        except ValueError:
            pass



def goterm2id(go_term_to_modif):
    return_id = int(go_term_to_modif.split(':')[1])
    return return_id


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


    import redis
    import gc
    from utils import config_utils




def _hog_lex_range(hog):
    """
    Decodes hog format
    :param hog: (bytes or string): hog id
    :return: hog_str: encoded hog id
    """
    hog_str = hog.decode() if isinstance(hog, bytes) else hog
    return hog_str.enco
def clearDB(dbnum):
    r = redis.StrictRedis(host='10.0.63.33', port=6379, db=dbnum)
    r.flushdb()
    return rde('ascii'), (hog_str[0:-1] + chr(1 + ord(hog_str[-1]))).encode('ascii')


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


def _get_go_terms(hog_id, oma, go_file):
    """
    Fetch the genes from hog id, then get all GO terms from it
    :param hog_id: hog id
    :return: list of GO term
    """
    genes = _get_hog_members(hog_id, oma)
    go_dict = {entry: {_format_go_term(e) for e in _get_entry_gene_ontology(oma, entry)} for entry in genes}
    go_dict_filtered = _filter_result(go_dict, go_file)

    return _clean_dictionary(go_dict_filtered)


def _format_go_term(e):
    # return e['TermNr']
    return 'GO:{:07d}'.format(e['TermNr'])


def _filter_result(go_dict, go_file):

    go_dict_filtered = {}

    for gene_name, terms in go_dict.items():
        filtered_terms = _filter_namespace(terms, go_file)
        if filtered_terms:
            go_dict_filtered[gene_name] = filtered_terms

    return go_dict_filtered


def _filter_namespace(list_terms, go_file, name='biological_process'):
    """
    Keep only go terms within the correct ontology
    :param list_terms: list of go terms
    :param name: namespace to keep; default: 'biological_process'
    :return: list of terms with the correct namespace
    """
    terms_to_keep = []
    for term in list_terms:
        try:
            if go_file[term].namespace == name:
                terms_to_keep.append(term)
        except KeyError:
            pass
    return terms_to_keep


def _clean_dictionary(dictionary):
    return {k: v for k, v in dictionary.items() if v}

if __name__ == '__main__':


    if preprocess_config.preprocessGO ==True:
        #Preprocess all of the GO terms' parents to avoid looking at the DAG
        obo_reader = obo_parser.GODag(obo_file=config_utils.datadir + 'project/data/go.obo')
        dt = h5py.special_dtype(vlen=np.dtype('int32'))
        omah5 = tables.open_file(config_utils.omadir + 'OmaServer.h5', mode='r')
        with h5py.File(config_utils.datadir + 'project/data/GOparents.h5', 'w', libver='latest') as h5_go_terms:
            start_time = time()
            h5_go_terms.create_dataset('goterms2parents', (10000000,), dtype=dt)
            dataset_go_terms_parents = h5_go_terms['goterms2parents']
            count = 0
            for go_term in obo_reader:
                go_term_read = obo_reader[go_term]
                if go_term_read.namespace == 'biological_process':

                    go_term_parents = go_term_read.get_all_parents()
                    go_term_parents_int = [goterm2id(go_term_read.id)] + [goterm2id(parent) for parent in go_term_parents]
                    dataset_go_terms_parents[goterm2id(go_term_read.id)] = go_term_parents_int
                    count += 1
                    if count % 1000 == 0:
                        print('saving')
                        h5_go_terms.flush()
            h5_go_terms.flush()
            print('Done with the parents in {} seconds'.format(time()-start_time))

    if preprocess_config.preprocessSTRINGDB:
        if preprocess_config.clearRedis == True:
            #clear the stringDB mapping
            StringRedisTOOLS.clearDB(1)

        r1 = redis.StrictRedis(host='localhost', port=6379, db=1)
        # sort the IDs alphanumerically.
        # protein1 protein2 neighborhood fusion cooccurence coexpression
        # experimental database textmining combined_score
        # 394.NGR_c00010 394.NGR_c33930 0 0 165 0 0 0 145 255
        # 37200000
        # save file line...
        refs = ['neighborhood', 'fusion', 'cooccurence', 'coexpression', 'experimental', 'database', 'textmining', 'combined_score']
        start_line = 0
        with open(preprocess_config.string_interactors +'/stringdata/protein.links.detailed.v10.5.txt', 'r') as stringAll:
            for i, line in enumerate(stringAll):
                if i > start_line:
                    words = line.split()
                    IDS = ''.join(sorted([words[0], words[1]]))
                    r1.set(IDS, stringAll.tell())
                if i % 1000000 == 0:
                    print(i)
    				if i % 10000000 == 0:
    					gc.collect()




    #annotate the HOG hash values h5 file with IDs from the uniprot mapper

    """
    to grab
    KEGG	KEGG_ID
    BioGrid	BIOGRID_ID
    ComplexPortal	COMPLEXPORTAL_ID
    DIP	DIP_ID
    STRING	STRING_ID

    """
    datasets = [ 'BIOGRID' , 'COMPLEXPORTAL' , 'DIP' , 'KEGG', 'STRING', 'OMA']
    #open uniprotmappings
    with open( config_utils.uniprotmappings , 'r')as uniprotmappings:
    #open hashes h5
        with h5py.File(config_utils.datadir +  , 'r+', libver='latest') as h5hashDB:
            for dataset in datasets:
                #add a datasets for each ID mapping in uniprot
                dt_2 = h5py.special_dtype(vlen=bytes)
                h5h5hashDB.create_dataset(dataset, (h5h5hashDB['hashes'].shape[0],) , dtype=dt_2)
            hog2goterms = h5_preprocess['hog2goterms']
            count = 1
            start_time = time()
            print('started!')
            oldID = ''
            mapdict ={}
            for i, row in enumerate(uniprotmappings):

                words= row.split()
                uniID = words[0]
                mapto = words[1]
                mapval = words[2]

                if mapto in datasets:
                    if mapto in mapdict:
                        mapdict[mapto].append(mapval)
                    else:
                        mapdic[mapto] =[ mapval]

                if uniID is not oldID:
                    if 'OMA' in mapdict:
                        entry = db_obj.entry_by_entry_nr(entry_nr)
                        fam = entry['OmaHOG']
                        for dataset in mapdict:
                            #add to annotation dictionary with each omaID
                            #to the correct family row in the hdf5 containing the hashes

                            if len(h5h5hashDB[dataset][fam]) > 0:
                                try:
                                    mapping = json.loads(h5h5hashDB[dataset][fam])
                                except ValueError:
                                    mapping = {}
                            else:
                                mapping = {}
                            mapping[mapdict['OMA']]=mapping[dataset]
                            h5h5hashDB[dataset][fam] = json.dumps(mapping)
                    mapdict = {}
                oldID = uniID

    #Add go terms from OMA
    with h5py.File(config_utils.datadir +  , 'r+', libver='latest') as h5hashDB:
        for row in h5h5hashDB['hashes']:
            #check if hash was compiled
            if len(row) > 0:
                hog_dict = _get_go_terms(fam2hogid(fam), omah5, obo_reader)
                hog2goterms[fam] = json.dumps(hog_dict).encode()
                if i % 1000 == 0:
                    print('saving {} {}'.format(time()-start_time, fam))
                    h5_go_terms.flush()


        #parse uniprot annotation File
        #each time you encounter an omaID get the hog

        #ADD the string info to each HOG
        #ADD the Kegg info to each HOG


        h5_go_terms.flush()


    print('DONE!')
