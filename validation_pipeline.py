import h5py

import ujson as json


def validate(hog_id):
    print(hog_id)


# load go terms data
parents_hdf5 = h5py.File('data/parents.h5', 'r')
hogs2goterms = parents_hdf5['hog2goterms']

# get hog with annotations
hogs_with_annotations = []
for fam, goterms in enumerate(hogs2goterms):
    if goterms.decode():
        if json.loads(goterms.decode()):
            hogs_with_annotations.append(fam)

            done_list = set()

# load profiler & LSH
# TODO


# for each hog with annotations, query results=
for hog in hogs_with_annotations:
    if hog not in done_list:
        ...
        do

        here call profiler

        stuff
        ...
        results_list = []
        for r in results_list:
            done_list.update(r)
        done_list.update(hog)


# display results