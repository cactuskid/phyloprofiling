from formatOrtho import convert_orthoxml_ids
import pyham

#pyham takes a file rather than a string, so save it as a local file
def retham(fam, dbObj, species_tree, datadir, replacement_dic, l=None):
    
    #dbObj.get_orthoxml(fam)
    #make ham analysis object
    #with tempfile.NamedTemporaryFile(dir = datadir ) as temp:
        #l.acquire()
    ortho = dbObj.get_orthoxml(fam)
    ortho = convert_orthoxml_ids(ortho, replacement_dic)
    hamObj = pyham.Ham(species_tree, ortho, type_hog_file="string")

    tree = pyhamtoTree(hamObj)

    return tree

        #l.release()
    #temp.write( ortho )

       # forget about temp file
       # return string of orthoxml
       # fam, dbobj, species tree, replacement dic
       # give obj to pyham with species tree; replace codes which do not work
       # pyham gibes back a ham object
       # coompute evo events (func in pyham)
       # returns pyham tree (ete3)

def yeildfams():
    for row in h5file.root.OrthoXML.Index:
        yield row[0]

def pyhamtoTree(hamOBJ):
	#use pyham to get all evolutionary events from a pyham object
	#turn into an array and a hash
	tp = pyham.TreeProfile(hamOBJ)
	tree = tp.compute_tree_profile_full()
	return tree