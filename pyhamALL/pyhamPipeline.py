from formatOrtho import convert_orthoxml_ids
import pyham

#pyham takes a file rather than a string, so save it as a local file
def get_ham(fam, dbObj, species_tree, replacement_dic, datadir=None, l=None):
	'''
	Get Ham object from fam id

	Args:
		fam : hog family
		dbObj : database object; read h5file
		species_tree : corrected species tree
		replacement_dic : mapping dictonary to correct orthoxml
		datadir = ??
		l = lock
	Returns :
		hamObj : pyham object
	'''
	

	#make ham analysis object
	#with tempfile.NamedTemporaryFile(dir = datadir ) as temp:
		#l.acquire()

	# get orthoXML from database
	ortho = dbObj.get_orthoxml(fam).decode()
	# correct orthoXML with mapping dict
	ortho = convert_orthoxml_ids(ortho, replacement_dic)
	# get ham Object

	hamObj = pyham.Ham(species_tree, ortho, type_hog_file="string", use_internal_name = True)

	return hamObj

	#l.release()
	#temp.write( ortho )


def yieldFamilies(h5file):
	'''
	Given a h5file containing OMA server, returns an iterator over the families
	'''
	for row in h5file.root.OrthoXML.Index:
		yield row[0]

def getOneFamily(i, h5file):
	'''
	get one family from database
	Args:
		i : family number
		h5file : OMA server file
	Return :
		family
	'''
	return h5file.root.OrthoXML.Index[i][0]

def pyhamtoTree(hamOBJ, fam):
	'''
	Use pyham to get all evolutionary events from a pyham object

	Args:
		hamOBJ : ham object, created from a species tree and an orthoXML file
		fam : hog family
	Returns :
		treemap : tree profile done from the ham object and the hog id (fam)
	'''
	hog = hamOBJ.get_hog_by_id(fam)
	tp = hamOBJ.create_tree_profile(hog=hog)

	return tp.treemap