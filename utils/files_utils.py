import ete3
import pandas as pd
from Bio import Entrez
#import config_utils


def get_tree(oma=None , overwrite = True ):
    if overwrite == True:
        ncbi = ete3.NCBITaxa()
        genome_ids_list = pd.DataFrame(oma.root.Genome.read())["NCBITaxonId"].tolist()
        tree = ncbi.get_topology(genome_ids_list)
        print(len(genome_ids_list))
        orphans = list(set(genome_ids_list) - set([int(x.name) for x in tree.get_leaves()]))
        print(len(orphans))
        Entrez.email = "clement.train@gmail.com"

        orphans_info = {}

        for x in orphans:
            search_handle = Entrez.efetch('taxonomy', id=str(x), retmode='xml')
            record = next(Entrez.parse(search_handle))
            orphans_info[x] = [x['TaxId'] for x in record['LineageEx']]

    tree = add_orphans(orphans_info, tree, genome_ids_list)

    for n in tree.traverse():
        if len( [ x for x in n.get_descendants()] ) == 1:
            # remove node with one Child
            parent = n.get_ancestors()[0]
            child = n.get_leaves()[0]
            parent.add_child(name=child.name)
            parent.add_child(name =n.name)
            child.delete()

    orphans = list(set(genome_ids_list) - set([int(x.name) for x in tree.get_leaves()]))
    print(orphans)
    tree_string = tree.write(format=1)

    return tree_string, tree


def generate_taxa_index(tree):
    """
    Generates an index for the global taxonomic tree for all OMA
    :param tree: ete3 tree
    :return: taxaIndex: dictionary key: node name (species name); value: index
        taxaIndexReverse: dictionary key: index: value: species name
    """
    taxa_index = {}
    taxa_index_reverse = {}
    for i, n in enumerate(tree.traverse()):
        taxa_index_reverse[i] = n.name
        taxa_index[n.name] = i

    return taxa_index, taxa_index_reverse


# def get_allowed_families(db_object, hog_level):
#
#     allowed_families_list = []
#     level = '''Level == b"{}"'''.format(hog_level)
#
#     for fam in db_object.get_hdf5_handle().get_node('/HogLevel').where(level):
#         allowed_families_list.append(fam[0])
#
#     return allowed_families_list


def add_orphans(orphan_info, tree, genome_ids_list, verbose=False):

    newdict = {}

    leaves = set([leaf.name for leaf in tree.get_leaves()])
    oldkeys = set(newdict.keys())
    orphans = set(genome_ids_list) - set([int(x.name) for x in tree.get_leaves()])
    print(orphans)
    keys = set()

    while len(orphans) > 0 and keys != oldkeys:
        oldkeys = keys
        for orphan in orphans:
            newdict[str(orphan_info[orphan][-1])] = str(orphan)
        keys = set(newdict.keys())

        for n in tree.traverse():
            if n.name in newdict and (newdict[n.name] not in leaves) and (n.name not in leaves):
                n.add_child(name=newdict[n.name])
                leaves.update(newdict[n.name])
                del newdict[n.name]

        orphans = set(genome_ids_list) - set([int(x.name) for x in tree.get_leaves()])
        print(orphans)
        for orphan in orphans:
            if len(orphan_info[orphan]) > 1:
                orphan_info[orphan].pop()
        newdict = {}

    return tree
