import ete3
import pandas as pd
from Bio import Entrez


def get_tree(oma=None , saveTree=True):
    ncbi = ete3.NCBITaxa()
    genome_ids_list = pd.DataFrame(oma.root.Genome.read())["NCBITaxonId"].tolist()
    genome_ids_list = [str(x) for x in genome_ids_list]
    tree = ncbi.get_topology(genome_ids_list , collapse_subspecies=False)

    print(len(genome_ids_list))
    orphans = list(set(genome_ids_list) - set([x.name for x in tree.get_leaves()]))
    print('missing taxa:')
    print(len(orphans))
    Entrez.email = "clement.train@gmail.com"
    orphans_info = {}
    for x in orphans:
        search_handle = Entrez.efetch('taxonomy', id=str(x), retmode='xml')
        record = next(Entrez.parse(search_handle))
        orphans_info[x] = [x['TaxId'] for x in record['LineageEx']]
    tree = add_orphans(orphans_info, tree, genome_ids_list)
    for n in tree.traverse():
        if len([x for x in n.get_descendants()]) == 1:
            # remove node with one Child
            parent = n.get_ancestors()[0]
            child = n.get_leaves()[0]
            parent.add_child(name=child.name)
            parent.add_child(name =n.name)
            child.delete()
    orphans = set(genome_ids_list) - set([x.name for x in tree.get_leaves()])
    tree_string = tree.write(format=1)
    if saveTree == True:
        with open( './mastertree.nwk' , 'w') as nwkout:
            nwkout.write(tree_string)

    print(orphans)

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


def add_orphans(orphan_info, tree, genome_ids_list, verbose=False):
    first = True


    newdict = {}

    leaves = set([leaf.name for leaf in tree.get_leaves()])

    orphans = set(genome_ids_list) - leaves
    oldkeys = set(list(newdict.keys()))

    keys = set()
    i = 0
    print(i)

    while first or ( len(orphans) > 0  and keys != oldkeys ) :
        first = False
        oldkeys = keys
        leaves = set([leaf.name for leaf in tree.get_leaves()])
        orphans = set(genome_ids_list) - leaves
        print(len(orphans))
        for orphan in orphans:
            if str(orphan_info[orphan][-1]) in newdict:
                newdict[str(orphan_info[orphan][-1])].append(orphan)
            else:
                newdict[str(orphan_info[orphan][-1])] = [orphan]
        keys = set(list(newdict.keys()))
        for n in tree.traverse():
            if n.name in newdict and n.name not in leaves:
                for orph in newdict[n.name]:
                    n.add_child(name=orph)
                del newdict[n.name]

        for orphan in orphans:
            if len(orphan_info[orphan]) > 1:
                orphan_info[orphan].pop()

        newdict = {}
    nodes = {}
    print(orphans)
    #clean up duplicates
    for n in tree.traverse():
        if n.name not in nodes:
            nodes[ n.name] =1
        else:
            nodes[ n.name] +=1

    for n in tree.traverse():
        if nodes[ n.name] >1:
            if n.is_leaf()== False:
                n.delete()
                nodes[ n.name]-= 1


    return tree
