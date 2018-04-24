import ete3
import re
import config
from tables import *
import pandas as pd

def create_species_tree(h5file, omaIdObj):
	'''
	Create and fix species tree; create a replacement dictionary used to remove special characters
	Args:
		h5file: OMA database
		omaIdObj: OMA id mapper from database object
	Returns:
		tree : species tree from NCBI
		replacement_dic : replacement dictionary used to correct the species tree and orthoxml files
	'''
	# load database
	# create ncbi object
	ncbi = ete3.NCBITaxa()
	# get genome list
	genome_ids_list = pd.DataFrame(h5file.root.Genome.read())["NCBITaxonId"].tolist()
	# get tree from NCBI; takes all the genomes ids and returns a species tree; some nodes are added between the given ids
	tree = ncbi.get_topology(genome_ids_list)
	# dictionary mapping NCBI taxa id with scientific names for all OMA genomes
	taxonId_SciName = {}
	for genome in h5file.root.Genome.read():
		taxonId_SciName[genome[0]] = genome[5].decode()
	
	# initialize replacement dictonary
	replacement_dic = {}
	# turns the tree into a string to parse it more easily
	tree_string = tree.write(format=1)
	# look for all species names with only 5 letters (special names from uniprot)
	uniprot_species = []
	uniprot_species = re.findall(r'\b[A-Z]{5}\b', tree_string)
	uniprot_species_to_add = ["STAA3","ECO57","BUCAI","CHLPN"]
	for species in uniprot_species_to_add:
		uniprot_species.append(species)
	# look for names from uniprot code for the special species names; store them in the replacement dictionary
	for species in uniprot_species:
		try:
			replacement_dic[species] = replace_characters(omaIdObj.genome_from_UniProtCode(species)[5].decode())
		except:
			pass
	# traverse the tree to fill the replacement dictionary and correct the 
	for node in tree.traverse():
		# the tree returned by ncbi contains more nodes than the provided genome id list, so the node as to be tested if its from the OMA database. 
		# If it is the case, the scientific name needs to be changed because OMA and NCBI use different notations
		if node.taxid in taxonId_SciName.keys():
			# replacement NCBI sci name by OMA name
			node.name = taxonId_SciName[node.taxid]
			# take the one from uniprot
			if node.name in uniprot_species:
				node.name = replacement_dic[node.name]
			# correct names; remove special characters breaking the species tree
			elif ',' in node.name or ('(' in node.name and ')'  in node.name) or ':' in node.name or '.' in node.name or ' ' in node.name:
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
	'''
	replace character from string
	Args:
		String
	Return:
		Corrected string
	'''
	for ch in ['.',',',' ','(',')',':']:
		if ch in string:
			string=string.replace(ch,'_')

	return string

def correct_orthoxml(instr, replacement_dic, verbose=False):
	'''Takes an orthoxml file as input, along with the replacement_dic, where the keys are scientific names (with 
	special characters already replaced) and values are the new name which matches the species tree.
	Replaces the special characters and old scientific names with the new ones.
	Args:
		instr: input string; orthoxml to correct
		replacement_dic: replacement dictionary used to correct orthoxml files and species tree
	Return:
		outstr: output string; corrected orthoxml
	'''
	outstr = ''
	exclude = []
	detected = True

	for line in instr.split('\n'):
		searchObj = re.search( r'.*<species name=\"(.*)\" NCBITaxId.*', line)

		if searchObj:
			old_name = searchObj.group(1)
			
			detected = False

			for key, value in replacement_dic.items():
				if old_name == key:
					line = line.replace(key, value)
					detected = True

			if verbose == True and detected == False:
				print(line)

		if detected == False:
			if '<gene id' in line:
				exclude += [s for s in line.split('"') if s.isdigit()]
				if verbose == True:
					print(exclude)

		if detected == True:
			if '<geneRef' in line:
				writeLine = True
				for ref in exclude:
					if ref in line:
						writeLine = False
				if writeLine == True:
					outstr += line + '\n' 
			else:
				outstr+= line + '\n'

	return outstr