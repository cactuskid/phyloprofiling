import pyham
import ete3
import xml.etree.ElementTree as ET


def get_orthoxml(fam, db_obj):
    orthoxml = db_obj.get_orthoxml(fam).decode()

    return orthoxml


def get_species_tree_from_orthoxml(orthoxml):
    species = get_species_from_orthoxml(orthoxml)
    print(len(species))
    ncbi = ete3.NCBITaxa()


    tree = ncbi.get_topology(species.keys())

    lineage = tree.get_linage(tree.get_tree_root().name)

    tree = ncbi.get_topology(lineage)

    # tree = ncbi.get_topology(ncbi.get_lineage(tree.name))


    print(tree)
    # prune tree, remove inter. node if only one child
    # correct nodes name
    for node in tree.traverse():
        # if len(node.children) == 1:
        #     print(node.name)
        #     node.delete()
        if node.name in species.keys():
            node.name = replacement_dic[species[node.name]]

    # turns the tree into a string

    tree_string = tree.write(format=1)


    return tree_string


def get_species_from_orthoxml(orthoxml):
    NCBITaxId2name = {}
    root = ET.fromstring(orthoxml)
    for child in root:
        if 'species' in child.tag:
            NCBITaxId2name[child.attrib['NCBITaxId']] = child.attrib['name']

    return NCBITaxId2name


# def get_ham_treemap_from_fam(fam, db_obj, species_tree, replacement_dic):
#     """
#     Get treemap ham object from fam id
#     :param fam: hog family
#     :param db_obj: database object; read h5file
#     :param species_tree: corrected species tree
#     :param replacement_dic: mapping dictionary to correct orthoxml
#     :return: tp.treemap
#     """
#     orthoxml = get_orthoxml(fam, db_obj, replacement_dic)
#     row = (fam, orthoxml)
#     species_tree = prune_tree(species_tree, orthoxml)
#     treemap = get_ham_treemap(row, species_tree)
#
#     return treemap


def get_ham_treemap_from_fam(fam, tree, db_obj):

    orthoxml = get_orthoxml(fam, db_obj)
    species_tree = get_species_tree_from_orthoxml(orthoxml, replacement_dic)

    ham_obj = pyham.Ham(species_tree, orthoxml.encode(), type_hog_file="orthoxml", use_internal_name=True,
                        orthoXML_as_string=True)
    hog = ham_obj.get_hog_by_id(fam)
    tp = ham_obj.create_tree_profile(hog=hog)

    return tp.treemap


def get_ham_treemap_from_row(row, tree, leaves):
    """
    Get treemap ham object from row (tuple: fam and orthoxml)
    :param row: tuple: fam and orthoxml
    :return: tp.treemap or none if fail
    """
    fam, orthoxml = row

    # TODO
    # check if tree is ok, repair it
    # TODO remove/change this
    species_tree = get_species_tree_from_orthoxml(orthoxml, replacement_dic)

    try:
        ham_obj = pyham.Ham(species_tree, orthoxml.encode(), type_hog_file="orthoxml",
                            use_internal_name=True, orthoXML_as_string=True)
        hog = ham_obj.get_hog_by_id(fam)
        tp = ham_obj.create_tree_profile(hog=hog)
        return tp.treemap
    except:
        return None


def yield_families(h5file, start_fam):
    """
    Given a h5file containing OMA server, returns an iterator over the families
    (not sure if still in use)
    :param h5file: omafile
    :param start_fam: fam to start on
    :return: fam number
    """
    for row in h5file.root.OrthoXML.Index:
        if row[0] > start_fam:
            yield row[0]


def get_one_family(i, h5file):
    '''
    get one family from database
    Args:
        i : family number
        h5file : OMA server file
    Return :
        family
    Not sure if still in use
    '''
    return h5file.root.OrthoXML.Index[i][0]


# def families2dataframe(h5file, verbose=False):
    # """
    # still in use ??
    # :param h5file:
    # :param verbose:
    # :return:
    # """
    # fams = h5file.root.OrthoXML.Index
    # df = dd.read_hdf(fams)
    # if verbose == True:
    #     print(df)
    # return df
