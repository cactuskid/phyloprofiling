
#orthoxml hack


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
			searchObj3 = re.search(r'.*<gene id=.*', line)
			if searchObj:
				old_name = searchObj.group(1)
				new_name = replace_characters(old_name)
				line = line.replace(old_name, new_name)

				for key, value in replacement_dic.items():
					if new_name == key:
						line = line.replace(key, value)

			if searchObj2:
				old_name = searchObj2.group(1)
				new_name = replace_characters(old_name)
				line = line.replace(old_name, new_name)

				for key, value in replacement_dic.items():
					if new_name == key:
						line = line.replace(key, value)
			if searchObj3:
				count = count + 1
			outstr+= line
	return outstr

def replace_characters(string):
	string = string.replace(".", "_")
	string = string.replace(" ", "_")
	string = string.replace("(", "")
	string = string.replace(")", "")
	string = string.replace(":", "-")
	return(string)

def fix_species_tree(species_tree):
	'''replaces characters which mess up the newick species tree'''
	replacement_dic = {}
	species_and_taxa = re.findall(r'"([^"]*)"', species_tree) 
	
	for old_name in species_and_taxa:
		new_name = replace_characters(old_name)
		replacement_dic[old_name] = new_name

	#sort by length of value so that long items get replaced first
	old_names_list = list(replacement_dic.keys())
	sorted_old_names_list = sorted(old_names_list , key = len, reverse=True)

	for name in sorted_old_names_list:
		species_tree = species_tree.replace( name, replacement_dic[name])

	species_tree = species_tree.replace("\"", "")
	species_tree = species_tree.replace("\n", "")
	
	return(species_tree)

def get_species_sciname(uniprotspeciescode):
	sciname = genome_df[genome_df['UniProtSpeciesCode']==uniprotspeciescode.encode("utf-8")]['SciName'].item().decode("utf-8") 
	return(sciname)

def generateHackedSpeciesTree(species_tree):
	species_tree = fix_species_tree(species_tree)
	with open( working_dir + 'speciestree_hack.nwk' , 'w') as outTree:
		outTree.write(species_tree)


#make a dictionary with the scientific name and the uniprot species code
#this replacement dic is for later to replace key with value in orthoxml
def generateReplacement(omaIdObj , species , verbose = False):



	t = Tree(species_tree, format =1 )
	print(t)
	species_tree = pyham.utils.get_newick_string(working_dir + "speciestree_hack.nwk", type="nwk")

	#get a list of all the species which are identified as their 5-letter uniprot species code in the species tree
	uniprot_species = []
	uniprot_species = re.findall(r'\b[A-Z]{5}\b', species_tree)
	uniprot_species.append("STAA3")
	uniprot_species.append("ECO57")
	uniprot_species.append("BUCAI")
	uniprot_species.append("CHLPN")

	print("UniProt 5-letter species codes which have replaced their scientific names in the species tree: "+str(uniprot_species))


	replacement_dic ={}
	for species in uniprot_species:
		try:
			replacement_dic[replace_characters(omaIdObj.genome_from_UniProtCode(species)[5].decode())] = species

		except:
			pass
	if verbose == True:
		print("\nreplacement_dic: "+ str(replacement_dic))
	return replacement_dic
