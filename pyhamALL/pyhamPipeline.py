import format_files
import pyham

#pyham takes a file rather than a string, so save it as a local file
def get_hamTree(fam, dbObj, species_tree, replacement_dic, l=None):
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
		treemap : tp.treemap
	'''
	#make ham analysis object
	if l:
		l.acquire()
	# get orthoXML from database
	ortho = dbObj.get_orthoxml(fam).decode()
	if l:
		l.realease()

	# correct orthoXML with mapping dict
	ortho = format_files.correct_orthoxml(ortho, replacement_dic, verbose=False)
	# get ham Object
	hamObj = pyham.Ham(species_tree, ortho.encode(), type_hog_file="orthoxml", use_internal_name = True, orthoXML_as_string=True)
	hog = hamObj.get_hog_by_id(fam)
	tp = hamObj.create_tree_profile(hog=hog)
	return tp.treemap


def yieldFamilies(h5file,startfam):
	'''
	Given a h5file containing OMA server, returns an iterator over the families
	'''
	for row in h5file.root.OrthoXML.Index:
		if row[0]>startfam:
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

def famsToDF(h5file, verbose=False):
	fams = h5file.root.OrthoXML.Index
	df = dd.read_hdf(fams)
	if verbose == True:
		print(df)
	return df
