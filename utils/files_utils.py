import ete3
import re
import pandas as pd

import pyham


def get_species_tree_replacement_dic(h5file, oma_id_obj):
    """
    Create and fix species tree; create a replacement dictionary used to remove special characters
    :param h5file: OMA database
    :param oma_id_obj: OMA id mapper from database object
    :return: tree: species tree from NCBI; replacement_dic: replacement dictionary used to correct the species tree and
    orthoxml files
    """
    # load database
    # create ncbi object
    ncbi = ete3.NCBITaxa()
    # get genome list
    genome_ids_list = pd.DataFrame(h5file.root.Genome.read())["NCBITaxonId"].tolist()
    # get tree from NCBI; takes all the genomes ids and
    # returns a species tree; some nodes are added between the given ids
    tree = ncbi.get_topology(genome_ids_list)
    # dictionary mapping NCBI taxa id with scientific names for all OMA genomes
    taxon_id_sci_name = {}
    for genome in h5file.root.Genome.read():
        taxon_id_sci_name[genome[0]] = genome[5].decode()

    # initialize replacement dictonary
    replacement_dic = {}
    # turns the tree into a string to parse it more easily
    tree_string = tree.write(format=1)
    # look for all species names with only 5 letters (special names from uniprot)
    uniprot_species = re.findall(r'\b[A-Z]{5}\b', tree_string)
    uniprot_species_to_add = ["STAA3", "ECO57", "BUCAI", "CHLPN"]
    for species in uniprot_species_to_add:
        uniprot_species.append(species)
    # look for names from uniprot code for the special species names; store them in the replacement dictionary
    for species in uniprot_species:
        try:
            replacement_dic[species] = replace_characters(oma_id_obj.genome_from_UniProtCode(species)[5].decode())
        except:
            pass
    # traverse the tree to fill the replacement dictionary and correct the
    for node in tree.traverse():
        # the tree returned by ncbi contains more nodes than the provided genome id list,
        # so the node as to be tested if its from the OMA database.
        # If it is the case, the scientific name needs to be changed because OMA and NCBI use different notations
        if node.taxid in taxon_id_sci_name.keys():
            # replacement NCBI sci name by OMA name
            node.name = taxon_id_sci_name[node.taxid]
            # take the one from uniprot
            if node.name in uniprot_species:
                node.name = replacement_dic[node.name]
            # correct names; remove special characters breaking the species tree
            elif ',' in node.name or \
                    ('(' in node.name and ')' in node.name) or ':' in node.name or '.' in node.name or ' ' in node.name:
                name = replace_characters(node.name)
                replacement_dic[node.name] = name
                node.name = name
        else:
            # the node is not present in OMA, keep NCBI notation
            node.name = node.sci_name

    # turns the tree into a string
    tree_fixed = tree.write(format=1)

    return replacement_dic, tree_fixed


def replace_characters(string):
    """
    Replace character from string
    :param string: string to correct
    :return: corrected string
    """
    for ch in ['.', ',', ' ', '(', ')', ':']:
        if ch in string:
            string = string.replace(ch, '_')

    return string


def correct_orthoxml(in_string, replacement_dic, verbose=False):
    """ Takes an orthoxml file as input, along with the replacement_dic, where the keys are scientific names (with
    special characters already replaced) and values are the new name which matches the species tree.
    Replaces the special characters and old scientific names with the new ones.
    :param in_string: input string; orthoxml to correct
    :param replacement_dic:
    :param verbose: replacement dictionary used to correct orthoxml files and species tree
    :return: output string; corrected orthoxml
    """
    output_string = ''
    exclude = []
    detected = True

    for line in in_string.split('\n'):
        search_object = re.search(r'.*<species name=\"(.*)\" NCBITaxId.*', line)

        if search_object:
            old_name = search_object.group(1)

            detected = False

            for key, value in replacement_dic.items():
                if old_name == key:
                    line = line.replace(key, value)
                    detected = True

            if verbose and not detected:
                print(line)

        if not detected:
            if '<gene id' in line:
                exclude += [s for s in line.split('"') if s.isdigit()]
                if verbose:
                    print(exclude)

        if detected:
            if '<geneRef' in line:
                write_line = True
                for ref in exclude:
                    if ref in line:
                        write_line = False
                if write_line:
                    output_string += line + '\n'
            else:
                output_string += line + '\n'

    return output_string


def get_ham_treemap(fam, db_obj, species_tree, replacement_dic):
    ortho = db_obj.get_orthoxml(fam).decode()

    ortho = correct_orthoxml(ortho, replacement_dic, verbose=False)
    ham_obj = pyham.Ham(species_tree, ortho.encode(), type_hog_file="orthoxml", use_internal_name=True,
                        orthoXML_as_string=True)
    hog = ham_obj.get_hog_by_id(fam)
    tp = ham_obj.create_tree_profile(hog=hog)

    return tp.treemap


def generate_taxa_index(species_tree):
    """
    Generates an index for the global taxonomic tree for all OMA
    :param species_tree: species tree in newick format
    :return: taxaIndex: dictionary key: node name (species name); value: index
        taxaIndexReverse: dictionary key: index: value: species name
    """
    t = ete3.Tree(species_tree, format=1)

    taxa_index = {}
    taxa_index_reverse = {}

    for i, node in enumerate(t.traverse()):
        taxa_index_reverse[i] = node.name
        taxa_index[node.name] = i
    return taxa_index, taxa_index_reverse


def get_allowed_families(db_object, hog_level):

    allowed_families_list = []
    level = '''Level == b"{}"'''.format(hog_level)

    for fam in db_object.get_hdf5_handle().get_node('/HogLevel').where(level):
        allowed_families_list.append(fam[0])

    return allowed_families_list
