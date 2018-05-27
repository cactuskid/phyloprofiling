from goatools import semantic


def resnik_sim_hdf5(go_id1, go_id2, godag, termcounts, hdf5):
    '''
        Computes Resnik's similarity measure.
    '''
    msca_goid = deepest_common_ancestor_hdf5([goterm2id(go_id1), goterm2id(go_id2)], godag, hdf5)
    return semantic.get_info_content(msca_goid, termcounts)


def deepest_common_ancestor_hdf5(go_ids, godag, hdf5):
    '''
        This function gets the nearest common ancestor
        using the above function.
        Only returns single most specific - assumes unique exists.
    '''
    # Take the element at maximum depth.
    return max(common_parent_go_ids_hdf5(go_ids, hdf5), key=lambda t: godag[t].depth)


def common_parent_go_ids_hdf5(go_ids, hdf5):
    '''
        This function finds the common ancestors in the GO
        tree of the list of goids in the input.
    '''

    candidates = set(hdf5['go_terms'][go_ids[0]].tolist())

    for go_id in go_ids[1:]:
        candidates_to_add = set(hdf5['go_terms'][go_id].tolist())

        candidates.intersection_update(candidates_to_add)

    corrected_candidates = [id2goterm(c) for c in candidates]
    return corrected_candidates


def goterm2id(go_term_to_modif):
    return int(go_term_to_modif.split(':')[1])


def id2goterm(go_term_to_modif):
    return 'GO:{:07d}'.format(go_term_to_modif)


