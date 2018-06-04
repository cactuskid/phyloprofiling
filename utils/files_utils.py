import ete3
import pandas as pd
from Bio import Entrez


def get_tree(oma):
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

    tree = addOrphans(orphans_info, tree, genome_ids_list)

    nodes = set([])


    corrected = True

    for n in tree.traverse():
        if len( [ x for x in n.get_descendants()] ) == 1:
            print(n.name)
            #remove node with one Child
            parent = n.get_ancestors()[0]
            print(parent)
            child = n.get_leaves()[0]
            parent.add_child(name = child.name)
            parent.add_child(name = n.name)
            child.delete()
            print(parent)
    count ={}
    for n in tree.get_leaves():
        if n.name not in count and n.name != '':
            count[n.name] = n

    tree.prune(list(count.values())+[tree])
    orphans = list(set(genome_ids_list) - set([int(x.name) for x in tree.traverse()]))

    print(orphans)

    tree = tree.write(format=1)

    return tree


def generate_taxa_index(h5file):
    """
    Generates an index for the global taxonomic tree for all OMA
    :param h5file: oma file
    :return: taxaIndex: dictionary key: node name (species name); value: index
        taxaIndexReverse: dictionary key: index: value: species name
    """
    taxa_index = {}
    taxa_index_reverse = {}
    for i, genome in enumerate(h5file.root.Genome.read()):
        taxa_index_reverse[i] = genome[5]
        taxa_index[genome[5]] = i

    return taxa_index, taxa_index_reverse


def get_allowed_families(db_object, hog_level):

    allowed_families_list = []
    level = '''Level == b"{}"'''.format(hog_level)

    for fam in db_object.get_hdf5_handle().get_node('/HogLevel').where(level):
        allowed_families_list.append(fam[0])

    return allowed_families_list


# def getParents(orphans, orthoxml,verbose):
#     # find stuff that is not in the species tree in the orthoxml and assign it to a taxonomic node
#
#     parentDict = {}
#     genes={}
#     root = ET.fromstring(orthoxml)
#     for elem in root:
#         if 'species' in elem.tag:
#             if elem.attrib['NCBITaxId'] in orphans:
#                 if elem.attrib['NCBITaxId'] not in genes:
#                     genes[elem.attrib['NCBITaxId']] =[]
#                 for gene in elem.iter():
#                     if 'gene' in gene.tag:
#                         try:
#                             genes[gene.attrib['id']] = elem.attrib['NCBITaxId']
#                         except KeyError:
#                             pass
#         if 'groups' in elem.tag:
#             parent_map = dict((c, p) for p in elem.getiterator() for c in p)
#             for groups in elem.iter():
#                 if 'geneRef' in groups.tag:
#                     if groups.attrib['id'] in genes:
#                         species = genes[groups.attrib['id']]
#                         if species not in parentDict:
#                             orthogroup = parent_map[groups]
#                             for prop in orthogroup:
#                                 if 'property' in prop.tag:
#                                     sciname = prop.get('value')
#                                     parentDict[sciname] = species
#                                     break
#
#     return parentDict


def addOrphans(orphan_info, tree, genome_ids_list, verbose=False):

    newdict = {}

    leaves = set([leaf.name for leaf in tree.get_leaves()])

    for orphan in orphan_info:
        newdict[str(orphan_info[orphan][-1])] = str(orphan)

    for n in tree.traverse():
        try:
            if n.name in newdict and newdict[n.name] not in leaves:
                n.add_child(name=newdict[n.name])
                leaves.add(newdict[n.name])
        except AttributeError:
            pass
    orphans = list(set(genome_ids_list) - set([int(x.name) for x in tree.get_leaves()]))
    oldkeys = set(newdict.keys())
    keys = ()

    while len(orphans) > 0 and keys != oldkeys:
        oldkeys = keys
        for orphan in orphans:
            if len(orphan_info[orphan])>1:
                orphan_info[orphan].pop()

            newdict[str(orphan_info[orphan][-1])] = str(orphan)
        keys = set(newdict.keys())

        for n in tree.traverse():
            try:
                if n.name in newdict and newdict[n.name] not in leaves:
                    n.add_child(name=newdict[n.name])
                    leaves.add(newdict[n.name])
            except AttributeError:
                pass

        orphans = list(set(genome_ids_list) - set([int(x.name) for x in tree.get_leaves()]))
    return tree


    # leftovers = set()
    #
    # if verbose:
    #     print(newdict)
    # for n in t.traverse():
    #     try:
    #         if n.sci_name in newdict and newdict[n.sci_name] not in leaves:
    #             n.add_child(name = newdict[n.sci_name])
    #             added.append(n.sci_name)
    #             leaves.add(newdict[n.sci_name])
    #     except AttributeError:
    #         pass
    #     # second attempt shortening the names...
    #     leftovers = set(newdict.keys()) - set(added)
    # if len(leftovers)>0:
    #     if verbose == True:
    #
    #         print('iterative start with leftovers:')
    #         print(leftovers)
    #
    #
    #     values = [ newdict[leftover] for leftover in leftovers]
    #     reduced = [ ''.join([word+' ' for word in leftover.split()[0:max(1,len(leftover.split())-1)]]).strip() for leftover in leftovers ]
    #     newdict = dict(zip(reduced,values))
    #
    #     reducedSet = set(reduced)
    #     reducedOld = set([])
    #
    #     while reducedSet != reducedOld :
    #         leaves = set([leaf.name for leaf in t.get_leaves()])
    #     if verbose:
    #         print('iterative start with leftovers:')
    #         print(leftovers)
    #     values = [newdict[leftover] for leftover in leftovers]
    #     reduced = [''.join([word+' ' for word in leftover.split()[0:max(1,len(leftover.split())-1)]]).strip() for leftover in leftovers ]
    #     newdict = dict(zip(reduced, values))
    #     reducedSet = set(reduced)
    #     reducedOld = set([])
    #     while reducedSet != reducedOld:
    #         for n in t.traverse():
    #             try:
    #                 if n.sci_name in newdict and newdict[n.sci_name] not in leaves:
    #                     n.add_child(name = newdict[n.sci_name])
    #                     leaves.add(newdict[n.sci_name])
    #                     added.append(n.sci_name)
    #                     if verbose:
    #                         print(n.sci_name)
    #             except AttributeError:
    #                 pass
    #
    #         leftoversNew = set(newdict.keys()) - set(added)
    #         if verbose:
    #             print(leftoversNew)
    #         if len(leftoversNew) ==0:
    #             if verbose:
    #                 print('DONE!')
    #             break
    #         values = [ newdict[leftover] for leftover in leftoversNew]
    #         reduced = [ ''.join([word+' ' for word in leftover.split()[0:max(1,len(leftover.split())-1)]]).strip() for leftover in leftoversNew]
    #         reducedOld = reducedSet
    #         reducedSet = set(reduced)
    #         newdict = dict(zip(reduced,values))
    #         if verbose:
    #             print('newdict')
    #             print(newdict)
    #
    #             print('newleftovers')
    #             print(leftoversNew)
    #
    # return t