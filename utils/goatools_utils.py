from goatools import semantic
import ujson as json
from utils import hashutils
from goatools.go_enrichment import GOEnrichmentStudy

##############enrichment##############################################

def return_enrichment_study_obj(gaf_taxfiltered):
    #make an enrichment study obj
    #obo_fname = download_go_basic_obo()
    obodag = GODag("go-basic.obo")
    goeaobj = GOEnrichmentStudy(
        gaf_taxfiltered.keys(), #
        gaf_taxfiltered, # geneid/GO associations possible with tree used for DB
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.15, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    return goeaobj

def run_GOEA_onresults(results, db_obj, goeaobj, outfile = None):
    #use lsh results to perform go enrichment
    #grab all ncbi ids
    geneids_study = [ member.omaid for member in [db_obj.iter_members_of_hog_id(int(result)) for result in results] ]
    goea_results_all = goeaobj.run_study(geneids_study)
    hogids =[ "HOG:" + (7-len(fam_id)) * '0' + fam_id for fam_id in HOGS[hog]['result'] ]
    prots = set([])
    for hogname in hogids:
        print(hogname)
        iterator = db_obj.iter_members_of_hog_id(hogname)
        omaids = frozenset([prot.omaid for prot in iterator ])
        HOGS[hog]['size']=len( omaids)
        HOGS[hog]['members'][hogname]=omaids
        prots = prots.union(omaids)
    goea_results_all = goeaobj.run_study(prots)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    goeaobj.wr_txt(folder + str(hog)+"enrichment.txt", goea_results_all)
    goea_results_terms = [ r.get_field_values(flds) for r in goea_results_sig ]
    goea_results_scores = [ r.p_fdr_bh for r in goea_results_sig ]
    HOGS[hog]['scores'] = goea_results_scores
    HOGS[hog]['terms'] = goea_results_terms
    HOGS[hog]['entrylist'] = list(prots)
    print(goea_results_terms)
    if outfile:
       with open(outfile , 'wb' ) as save:
           save.write(pickle.dumps(HOGS,2))




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
        Fins the common ancestors in the GO
        tree of the list of goids in the input.
    '''
    candidates = set(hdf5_set[go_ids[0]].tolist())
    for go_id in go_ids[1:]:
        candidates_to_add = set(hdf5_set[go_id].tolist())
        candidates.intersection_update(candidates_to_add)
    corrected_candidates = [id2goterm(c) for c in candidates]
    return corrected_candidates

def get_go_terms_hdf5(hog_id, hdf5_set):
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
