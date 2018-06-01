import pyham
import xml.etree.cElementTree as ET


def get_orthoxml(fam, db_obj):
    orthoxml = db_obj.get_orthoxml(fam).decode()

    return orthoxml


# def get_species_tree_from_orthoxml(orthoxml):
#     species = get_species_from_orthoxml(orthoxml)
#     print(len(species))
#     ncbi = ete3.NCBITaxa()
#
#
#     tree = ncbi.get_topology(species.keys())
#
#     lineage = tree.get_linage(tree.get_tree_root().name)
#
#     tree = ncbi.get_topology(lineage)
#
#     # tree = ncbi.get_topology(ncbi.get_lineage(tree.name))
#
#
#     print(tree)
#     # prune tree, remove inter. node if only one child
#     # correct nodes name
#     for node in tree.traverse():
#         # if len(node.children) == 1:
#         #     print(node.name)
#         #     node.delete()
#         if node.name in species.keys():
#             node.name = replacement_dic[species[node.name]]
#
#     # turns the tree into a string
#
#     tree_string = tree.write(format=1)
#
#
#     return tree_string


def get_species_from_orthoxml(orthoxml):
    NCBITaxId2name = {}
    root = ET.fromstring(orthoxml)
    for child in root:
        if 'species' in child.tag:
            NCBITaxId2name[child.attrib['NCBITaxId']] = child.attrib['name']

    return NCBITaxId2name

def switch_name_ncbiid(orthoxml):

    root = ET.fromstring(orthoxml)
    #orthoxml = orthoxml.decode('utf-8')
    for child in root:
        if 'species' in child.tag:
            orthoxml = orthoxml.replace(child.attrib['name'], child.attrib['NCBITaxId'])

    return orthoxml#.encode()


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


def get_species_tree_from_orthoxml(orthoxml , tree, leaves , verbose = False):
    # configure this function using a partial and give it the tree and the set of all leaf names
    # only adjust the tree when there is stuff in the orthoxml that isnt in the tree
    species = get_species_from_orthoxml(orthoxml)
    orphans = (set(species) - leaves)
    if len(orphans) == 0:
        return tree.write(format=1)
    else:
        parents = getParents(orphans, orthoxml , verbose)
        tree = addOrphans( parents , tree, verbose)
        orphans = (set(species)- set( [node.name for node in tree.get_leaves()]))
        if verbose:
            print(orphans)
        return tree.write(format=1)


def getParents(orphans, orthoxml , verbose):
    # find stuff that is not in the species tree in the orthoxml and assign it to a taxonomic node
    parentDict = {}
    genes={}
    root = ET.fromstring(orthoxml)
    for elem in root:
        if 'species' in elem.tag:
            if elem.attrib['NCBITaxId'] in orphans:
                if elem.attrib['NCBITaxId'] not in genes:
                    genes[elem.attrib['NCBITaxId']] =[]
                for gene in elem.iter():
                    if 'gene' in gene.tag:
                        try:
                            genes[gene.attrib['id']] = elem.attrib['NCBITaxId']
                        except KeyError:
                            pass
        if 'groups' in elem.tag:
            parent_map = dict((c, p) for p in elem.getiterator() for c in p)
            for groups in elem.iter():
                if 'geneRef' in groups.tag:
                    if groups.attrib['id'] in genes:
                        species = genes[groups.attrib['id']]
                        if species not in parentDict:
                            orthogroup = parent_map[groups]
                            for prop in orthogroup:
                                if 'property' in prop.tag:
                                    sciname = prop.get('value')
                                    parentDict[sciname] = species
                                    break
    return parentDict


def addOrphans(parentDict, t, verbose=False):
    # add orphans to tree
    added =[]
    newdict = parentDict
    leftovers = set()
    if verbose:
        print(newdict)
    for n in t.traverse():
        try:
            if n.sci_name in newdict:
                n.add_child(name = newdict[n.sci_name])
                added.append(n.sci_name)
        except AttributeError:
            pass
        # second attempt shortening the names...
        leftovers = set(newdict.keys()) - set(added)
    if len(leftovers)>0:
        if verbose:
            print('iterative start with leftovers:')
            print(leftovers)
        values = [newdict[leftover] for leftover in leftovers]
        reduced = [''.join([word+' ' for word in leftover.split()[0:max(1,len(leftover.split())-1)]]).strip() for leftover in leftovers ]
        newdict = dict(zip(reduced, values))
        reducedSet = set(reduced)
        reducedOld = set([])
        while reducedSet != reducedOld:
            for n in t.traverse():
                try:
                    if n.sci_name in newdict:
                        n.add_child(name = newdict[n.sci_name])
                        added.append(n.sci_name)
                        if verbose:
                            print(n.sci_name)
                except AttributeError:
                    pass
            leftoversNew = set(newdict.keys()) - set(added)
            if verbose:
                print(leftoversNew)
            if len(leftoversNew) ==0:
                if verbose:
                    print('DONE!')
                break
            values = [ newdict[leftover] for leftover in leftoversNew]
            reduced = [ ''.join([word+' ' for word in leftover.split()[0:max(1,len(leftover.split())-1)]]).strip() for leftover in leftoversNew]
            reducedOld = reducedSet
            reducedSet = set(reduced)
            newdict = dict(zip(reduced,values))
            if verbose:
                print('newdict')
                print(newdict)
                print('newleftovers')
                print(leftoversNew)
    return t


def get_ham_treemap_from_fam(fam, tree, db_obj):
    orthoxml = get_orthoxml(fam, db_obj)
    orthoxml = switch_name_ncbiid(orthoxml)
    row = (fam, orthoxml)
    ham_obj = pyham.Ham(tree, orthoxml.encode(), type_hog_file="orthoxml", use_internal_name=False,
                        orthoXML_as_string=True)
    hog = ham_obj.get_hog_by_id(fam)
    tp = ham_obj.create_tree_profile(hog=hog)
    return tp.treemap


def get_ham_treemap_from_row(row, tree , leaves):
    fam, orthoxml = row
    treestr = get_species_tree_from_orthoxml(orthoxml, tree, leaves , verbose = False)
    # try:
    ham_obj = pyham.Ham(treestr, orthoxml, type_hog_file="orthoxml", use_internal_name=True, orthoXML_as_string=True)
    hog = ham_obj.get_hog_by_id(fam)
    tp = ham_obj.create_tree_profile(hog=hog)
    return tp.treemap
    # except:
    #     return None


# def get_ham_treemap_from_row(row, tree, leaves):
#     """
#     Get treemap ham object from row (tuple: fam and orthoxml)
#     :param row: tuple: fam and orthoxml
#     :return: tp.treemap or none if fail
#     """
#     fam, orthoxml = row
#
#     # TODO
#     # check if tree is ok, repair it
#     # TODO remove/change this
#     species_tree = get_species_tree_from_orthoxml(orthoxml, replacement_dic)
#
#     try:
#         ham_obj = pyham.Ham(species_tree, orthoxml.encode(), type_hog_file="orthoxml",
#                             use_internal_name=True, orthoXML_as_string=True)
#         hog = ham_obj.get_hog_by_id(fam)
#         tp = ham_obj.create_tree_profile(hog=hog)
#         return tp.treemap
#     except:
#         return None


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
