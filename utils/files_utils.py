import ete3
import pandas as pd


def get_tree(oma):
    ncbi = ete3.NCBITaxa()
    genome_ids_list = pd.DataFrame(oma.root.Genome.read())["NCBITaxonId"].tolist()
    tree = ncbi.get_topology(genome_ids_list)
	for n in tree.traverse():
		if n.get_descendants() ==1:
			#remove node with one Child
			n.delete()
    return tree


def get_leaves(newick_tree):
    leaves = set(newick_tree.get_leaves())

    return leaves


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
