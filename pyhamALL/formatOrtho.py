import ete3
import re
import config
import pyham
#orthoxml hack
import fnmatch

def convert_orthoxml_ids(instr, replacement_dic):
	'''Takes an orthoxml file as input, along with the replacement_dic, where the keys are scientific names (with 
	special characters already replaced) and values are the new name which matches the species tree.
	Replaces the special characters and old scientific names with the new ones.
	Returns the number of genes in the orthoxml file.'''
	mylist = []
	count = 0 
	outstr = ''
	for line in instr.split('\n'):
		searchObj = re.search( r'.*<species name=\"(.*)\" NCBITaxId.*', line)
		searchObj2 = re.search(r'.*<property name=\"TaxRange\" value=\"(.*)\"\/>', line)
		if searchObj:
			old_name = searchObj.group(1)
			#new_name = replace_characters(old_name)
			#line = line.replace(old_name, new_name)

			for key, value in replacement_dic.items():
				if old_name == key:
					line = line.replace(key, value)

		if searchObj2:
			old_name = searchObj2.group(1)
			#new_name = replace_characters(old_name)
			#line = line.replace(old_name, new_name)

			for key, value in replacement_dic.items():
				if old_name == key:
					line = line.replace(key, value)
		
		outstr+= line
	return outstr

def replace_characters(string):
	'''
	replace character from string
	Args:
		String
	Return:
		Corrected string
	'''
	for ch in ['.',' ','(',')',':']:
		if ch in string:
			string=string.replace(ch,'_')

	return string



def get_species_sciname(uniprotspeciescode):
	sciname = genome_df[genome_df['UniProtSpeciesCode']==uniprotspeciescode.encode("utf-8")]['SciName'].item().decode("utf-8") 
	return(sciname)

def generateHackedSpeciesTree(species_tree):
	species_tree = fix_species_tree(species_tree)
	with open(config.working_dir + 'speciestree_hack.nwk' , 'w') as outTree:
		outTree.write(species_tree)

	return species_tree

#make a dictionary with the scientific name and the uniprot species code
#this replacement dic is for later to replace key with value in orthoxml
def fix_species_tree(species_tree, omaIdObj , verbose = False):
	'''replaces characters which mess up the newick species tree

		Args:
			species_tree : species tree in newick format
			omaIdObj : Oma ID object
			verbose : boolean,  output more info

	'''

	# need two replacement dictionary
	# one is returned and one is used to correct the tree

	replacement_dict_species_tree = {}
	replacement_dict_orthoXML = {}
	
	#species_and_taxa = re.findall(r'"([^"]*)"', species_tree) 
	
	t = ete3.Tree( species_tree, format=1 , quoted_node_names = True )
	
	for node in t.traverse():
		if verbose == True:
			if node.is_root():
				print('name')
				print(node.name)


		if ',' in node.name or ('(' in node.name and ')'  in node.name) or ':' in node.name or '.' in node.name:
			name = replace_characters(node.name)
			replacement_dict_orthoXML[node.name] = name 
			

	t = ete3.Tree( t.write(format=1) , format=1)
	# create a second tree, is a copy from the first one, but has LUCA for root
	t2 = ete3.Tree()
	t2.name="LUCA"
	for n in t.get_children():
		t2.add_child(n)
	for node in t2.traverse():
		if verbose == True:
			if node.is_root():
				print('root name')
				print(node.name)
				# create second replacement dictionary
		if ',' in node.name or ('(' in node.name and ')'  in node.name) or ':' in node.name or '.' in node.name:
			name = replace_characters(node.name)
			replacement_dict_species_tree[node.name] = name 		
	species_tree = t2.write(format=1, format_root_node= True)
	
	if verbose == True:
		print(species_tree.find('LUCA'))
	
	#species_tree = species_tree.replace("\n", "")
	
	#get a list of all the species which are identified as their 5-letter uniprot species code in the species tree
	uniprot_species = []
	uniprot_species = re.findall(r'\b[A-Z]{5}\b', species_tree)

	uniprot_species_to_add = ["STAA3","ECO57","BUCAI","CHLPN"]
	for species in uniprot_species_to_add:
		uniprot_species.append(species)

	if verbose == True:
		print("UniProt 5-letter species codes which have replaced their scientific names in the species tree: "+str(uniprot_species))
	
	for species in uniprot_species:

		try:
			replacement_dict_species_tree[species] = replace_characters(omaIdObj.genome_from_UniProtCode(species)[5].decode())
			replacement_dict_orthoXML[species] = replace_characters(omaIdObj.genome_from_UniProtCode(species)[5].decode())
		except:
			pass
		
	if verbose == True:
		#print("\nreplacement_dic: "+ str(replacement_dic))
		for species in uniprot_species:
			if species in replacement_dict_species_tree:
				print( replacement_dict_species_tree[species]) 
	
	for name in replacement_dict_species_tree:
		species_tree = species_tree.replace(name, replacement_dict_species_tree[name])

	species_tree = species_tree.replace("\"", "")
	
	with open(config.working_dir + 'speciestree_hack.nwk' , 'w') as outhack:
		outhack.write(species_tree)


	return species_tree , replacement_dict_orthoXML
