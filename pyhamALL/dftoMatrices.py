#use dask and pyham to create a big sparse matrix for each type of evolutionary event
import pyham
import dask
import datasketch


def pyhamtoArray(hamOBJ):
	#use pyham to get all evolutionary events from a pyham object
	#turn into an array and a hash
	
	