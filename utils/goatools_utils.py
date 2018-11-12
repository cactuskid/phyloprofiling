from goatools import semantic
from goatools.obo_parser import GODag

import ujson as json
from utils import hashutils
from utils import config_utils
import pickle
from goatools.go_enrichment import GOEnrichmentStudy

##############enrichment##############################################

def return_enrichment_study_obj(gaf_taxfiltered):
    '''
    Generate go enrichment study object with a background dataset.
    '''

    obodag = GODag(config_utils.datadir+"/GOData/go-basic.obo")
    goeaobj = GOEnrichmentStudy(
        gaf_taxfiltered.keys(), #
        gaf_taxfiltered, # geneid/GO associations possible with tree used for DB
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.15, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    return goeaobj

def buildGAF(gaf_file , universe= None):
    gaf_filtered = {}
    with open(gaf_file, mode='r') as gafin:
        for line in gafin:
            words = line.split()
            if words[0] not in gaf_filtered:
                gaf_filtered[words[0]]=set([words[1]])
            else:
                gaf_filtered[words[0]].add(words[1])

    if universe:
        gaf_filtered = { prot:gaf_filtered[prot] for prot in universe}


    return gaf_filtered

def run_GOEA_onresults(results, db_obj, goeaobj, outname = None):
    '''
        Perform enrichment analysis on returned results
        grabs all member protein of all hogs in result
        returns goe results and HOG composition
    '''

    #print(db_obj.member_of_hog_id(int(results[0])))
    hogids =[ "HOG:" + (7-len(fam_id)) * '0' + fam_id for fam_id in results ]
    #print( db_obj.member_of_hog_id(hogids[0]) )
    HOGS={}
    print('compiling hogs')
    prots = []

    for i,result in enumerate(hogids):
        if i %10 ==0:
            print(i)

        HOGS[result]=[]

        for member in db_obj.iter_members_of_hog_id(result):
            HOGS[result].append(member.omaid)
            prots.append(member.omaid)
    print('done')
    print('running GO enrichment study')


    goea_results_all = goeaobj.run_study(prots )
    print('done')
    with open( config_utils.datadir + outname + 'Hogs2Prots.pkl' , 'wb' ) as save:
       save.write(pickle.dumps(HOGS,2))

    goeaobj.wr_txt(config_utils.datadir+ str(outname)+"enrichment.txt", goea_results_all)
    print('DONE!')
    return goea_results_all, HOGS






######################resnik semsim ###################################################

def resnik_sim_hdf5(go_id1, go_id2, godag, termcounts, hdf5):
    '''
        Computes Resnik's similarity measure.
    '''
    try:
        msca_goid = deepest_common_ancestor_hdf5([goterm2id(go_id1), goterm2id(go_id2)], godag, hdf5)
        score = semantic.get_info_content(msca_goid, termcounts)
    except:
        score = -1
    return score

def deepest_common_ancestor_hdf5(go_ids, godag, hdf5):
    '''
        Gets the nearest common ancestor
        using the above function.
        Only returns single most specific - assumes unique exists.
    '''
    # Take the element at maximum depth.
    return max(common_parent_go_ids_hdf5(go_ids, hdf5), key=lambda t: godag[t].depth)


def common_parent_go_ids_hdf5(go_ids, hdf5_set):
    '''
        Finds the common ancestors in the GO
        tree of the list of goids in the input.
    '''
    candidates = set(hdf5_set[go_ids[0]].tolist())
    for go_id in go_ids[1:]:
        candidates_to_add = set(hdf5_set[go_id].tolist())
        candidates.intersection_update(candidates_to_add)
    corrected_candidates = [id2goterm(c) for c in candidates]
    return corrected_candidates

def get_go_terms_hdf5(hog_id, hdf5_set):
    '''
        grabs the preprocessed GO term info for HOGs from an HDF5.
    '''
    fam = hashutils.hogid2fam(hog_id)
    try:
        go_terms = json.loads(hdf5_set[fam])
    except:
        go_terms = None
    return go_terms

def goterm2id(go_term_to_modif):

    return int(go_term_to_modif.split(':')[1])

def id2goterm(go_term_to_modif):
    return 'GO:{:07d}'.format(go_term_to_modif)
