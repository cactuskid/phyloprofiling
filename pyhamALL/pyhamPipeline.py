#pyham takes a file rather than a string, so save it as a local file
def retham(fam, l, dbObj, species_tree, datadir , replacement_dic):
    dbObj.get_orthoxml(fam)
    #make ham analysis object
    with tempfile.NamedTemporaryFile(dir = datadir ) as temp:
        l.acquire()
        ortho = dbObj.get_orthoxml(fam)
        l.release()
        temp.write( ortho )

def yeildfams():
    for row in h5file.root.OrthoXML.Index:
        yield row[0]