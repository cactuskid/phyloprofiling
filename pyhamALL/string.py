#string db api

import json
import urllib as urllib
import networkx as nx
import xml.etree.ElementTree as ET

def call_string_list(species, identifiers, required_score,additional_network_nodes, verbose = False):
	
	#make a networkx graph from some IDs 
	#this functions reurns the xml and a list of seeds

	adress = 'http://string-db.org/api/'
	resolve='tsv-no-header/resolve?'
	
	#construct args
	#map Ids to stringDB
	identifier_str ='identifier='
	format_str ='format=only-ids'
	seeds ={}
	for prot in identifiers:		
		args = [ identifier_str + prot , format_str]
		url = adress + resolve
		for arg in args:
			if arg != args[-1]:
				url = url + arg + '&'
			else:
				url = url + arg
		
		seeds[prot] = []
		#get stringdb ids
		ids = urllib.urlopen(url).read()
		ids = ids.split('\n')
		for uniprot in ids:
			if uniprot != '' and 'not found' not in uniprot and 'Error' not in uniprot and 'Generally' not in uniprot:
				if prot in seeds.keys():
					seeds[prot].append(uniprot)
				else:
					seeds[prot]=[]
					seeds[prot].append(uniprot)
	graphs = {}

	#once youve mapped then pull down the net

	if len(seeds)>0:
		#get ineractions and interactors
		queryInteraction = 'psi-mi/interactions?'
		queryInteractor = 'psi-mi/interactors?'
		score_str = 'required_score='+ str(required_score)
		additional_str= 'additional_network_nodes='+str(additional_network_nodes)
		identifier_str ='identifier='
		
		for seed in seeds.keys():
			for prot in seeds[seed]:
				args = [ identifier_str + prot, score_str, additional_str]
				url = adress + queryInteraction
				for arg in args:
					if arg != args[-1]:
						url = url + arg + '&'
					else:
						url = url + arg
				if verbose == True:
					print( 'loading interactions for ' +prot)
					print(url)
				try:
				 	graphs[seed]
				except:
					graphs[seed] = []
				graphs[seed].append( urllib.urlopen(url).read())
	return graphs
	
def stringtograph(xmlfile, seeds , G= None , verbose = False):
	#you can iteratively build a graph by just refeeding it to this function with new data
	
	if G is None:
		G=nx.Graph()

	try:
		tree = ET.parse(xmlfile)
		root = tree.getroot()
		element = root.find('{net:sf:psidev:mi}entry')
		interactors = element.find('{net:sf:psidev:mi}interactorList')
		interactions = element.find('{net:sf:psidev:mi}interactionList')

	except:
		print('xml parse error')
	#keep uniport refs
	dic={}
	try:
		for child in interactors.iter('{net:sf:psidev:mi}interactor'):
			
				ID = child.get('id')
				xref = child.find('{net:sf:psidev:mi}xref')
				primaryRef = xref.find('{net:sf:psidev:mi}primaryRef')
				name = primaryRef.get('id')
				references = {}
				refiter = xref.findall('{net:sf:psidev:mi}secondaryRef')
				for ref in refiter:
					if ref.get('db') not in references.keys():
						references[ref.get('db')] = []
					references[ref.get('db')].append(ref.get('id'))
			
				if verbose ==True:
					print( name)
				if name in seeds:
					notSeed= False
					try:
						G.node[name]
						x = G.node[name]['ref'].copy()
						x.update(references)
						G.node[name]['ref'] = x
						dic[ID] = name
					except:
						G.add_node(name ,id =name)
						G.node[name]['ref'] = references
						G.node[name]['seed'] = True
						G.node[name]['color'] = 'red'
						dic[ID] = seedref
				else:
					G.add_node(name, id = name)
					G.node[name]['ref'] = references
					G.node[name]['seed'] = False
					G.node[name]['color'] = 'blue'
					dic[ID] = name
	except:
			print( 'xml element tree error: couldnt find interactors')
			return G


	for child in interactions:
		#this part is sloppier than the rest... 
		#gonna add filters for which evidence is used
		
		try:
			ID = child.attrib['id']
		except:
			print( 'xml element tree error: couldnt find interactions')
		
		participants = []
		try:
			for participant in child.find('{net:sf:psidev:mi}participantList'):
				node = participant.find('{net:sf:psidev:mi}interactorRef').text
				participants.append(dic[node])
			if len(participants)>1:
				G.add_edge(participants[0], participants[1])
				
				#todo add experiment data to each link

				references = []
				explist = child.find('{net:sf:psidev:mi}experimentList')
				for ref in explist.iter():
					if ref.text not in references:
						references.append(ref.text)
				G[participants[0]][participants[1]]['ref'] = references
		except:
			print('couldnt find')
	if verbose ==True:
		print( str(len(G.nodes())) + ' nodes')
		print( str(len(G.edges())) + ' edges')
	return G


		

	

