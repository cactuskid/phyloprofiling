from tables import *
import functools
import pandas as pd
import multiprocessing as mp
import time
import h5sparse
import pickle
from datasketch import MinHashLSH, MinHashLSHForest
import scipy

from pyoma.browser import db

from utils import config_utils
from utils import files_utils
from utils import pyhamutils


def generate_dataframes(size=100):
    families = {}
    with open_file(config_utils.omadir + 'OmaServer.h5', mode="r") as h5OMA:
        db_obj = db.Database(h5OMA)
        oma_id_obj = db.OmaIdMapper(db_obj)
        dic, tree = files_utils.get_species_tree_replacement_dic(h5OMA, oma_id_obj)
        # taxaIndex, reverse = files_utils.generate_taxa_index(tree)
        READ_ORTHO = functools.partial(pyhamutils.get_orthoxml, dbObj=db_obj, species_tree=tree, replacement_dic=dic)
        for i, row in enumerate(h5OMA.root.OrthoXML.Index):
            fam = row[0]
            families[fam] = {'ortho': READ_ORTHO(fam)}
            if len(families) > size:
                pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                pd_dataframe['Fam'] = pd_dataframe.index
                families = {}
                yield pd_dataframe


def worker(i, q, retq, matq, l):
    print('worker init ' + str(i))
    while True:
        df = q.get()
        if df is not None:



            df['tree'] = df[['Fam', 'ortho']].apply(HAMPIPELINE, axis=1)
            df['hash'] = df[['Fam', 'tree']].apply(HASHPIPEline, axis=1)

            # df['rows'] = df['tree'].apply( ROWPIPELINE )
            retq.put(df[['Fam', 'hash']])
            # matq.put(df['rows'])
        else:
            print('DONE WORKER' + str(i))
            break


def matrixupdater(i, q, retq, matq, l, rows, columns):
    hogmat = scipy.sparseCSR((rows, columns))
    printstart = time.clock()
    savestart = time.clock()
    globaltime = time.clock()

    datestr = "{:%B_%d_%Y_%H_%M}".format(dt.now())

    while True:
        rows = matq.get()
        if rows != None:
            for fam, sparserow in rows.itesm():
                hogmat[fam, :] = sparserow
            if time.clock() - savestart > 2000:
                with h5sparse.File(config_utils.datadir + datestr + "matrix.h5", 'w') as h5matrix:
                    h5matrix.create_dataset('hogmat', data=hogmat)

        else:

            with h5sparse.File(config_utils.datadir + datestr + "matrix.h5", 'w') as h5matrix:
                h5matrix.create_dataset('hogmat', data=hogmat)
            print('DONE MAT UPDATER' + str(i))

            break


def saver(i, q, retq, matq, l):
    dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']

    printstart = time.clock()
    savestart = time.clock()
    globaltime = time.clock()

    chunksize = 100
    count = 0
    dataset_names = ['Fam', 'duplication', 'gain', 'loss', 'presence']
    threshold = 0.7
    lsh = MinHashLSH(threshold=threshold, num_perm=128)
    forest = MinHashLSHForest(num_perm=128)

    dt = datetime
    datestr = "{:%B_%d_%Y_%H_%M}".format(dt.now())

    with open('errors.txt', 'w')as errfile:

        with h5py.File(config_utils.datadir + datestr + 'hashes.h5', 'w', libver='latest') as h5hashes:
            dsets = {}
            for dataset_name in dataset_names:
                print(dataset_name)
                if dataset_name not in list(h5hashes.keys()):
                    dataset = h5hashes.create_dataset(dataset_name, (chunksize, 0), maxshape=(None, None),
                                                      dtype='int32')
                dsets[dataset_name] = h5hashes[dataset_name]
            print(dsets)

            print('saver init ' + str(i))
            while True:
                thisdf = retq.get()
                print(time.clock() - globaltime)
                print(count)
                if thisdf is not None:

                    hashes = thisdf['hash'].to_dict()
                    for fam in hashes:
                        if hashes[fam] is not None:

                            for famhashname in hashes[fam]['dict']:
                                lsh.insert(famhashname, hashes[fam]['dict'][famhashname])
                                forest.add(famhashname, hashes[fam]['dict'][famhashname])
                            hashvals = hashes[fam]['hashes']

                            for i, event in enumerate(hashvals):
                                dataset = dataset_names[i]
                                if len(dsets[dataset]) < fam + 10:
                                    dsets[dataset].resize((fam + chunksize, len(hashvals[event])))
                                dsets[dataset][fam, :] = hashvals[event]
                            h5hashes.flush()
                        else:
                            print('error')
                            print(fam)
                            errfile.write(str(fam) + '\n')

                    if time.clock() - printstart > 60:
                        print(thisdf['Fam'].max())
                        print(time.clock() - globaltime)
                        printstart = time.clock()

                    if time.clock() - savestart > 2000:
                        print(thisdf['Fam'].max())
                        print('saving')
                        with open(config_utils.datadir + datestr + '_' + str(threshold) + '_' + 'newlsh.pkl', 'wb') as lshout:
                            with open(config_utils.datadir + datestr + 'newlshforest.pkl', 'wb') as forestout:
                                pickle.dump(lsh, lshout, -1)
                                pickle.dump(forest, forestout, -1)
                        savestart = time.clock()
                        print(time.clock() - globaltime)
                    count += len(thisdf)
                else:
                    with open(config_utils.datadir + datestr + '_' + str(threshold) + '_' + 'newlsh.pkl', 'wb') as lshout:
                        with open(config_utils.datadir + datestr + 'newlshforest.pkl', 'wb') as forestout:
                            pickle.dump(lsh, lshout, -1)
                            pickle.dump(forest.index(), forestout, -1)

                    print('DONE UPDATER' + str(i))
                    break


def mp_with_timeout(number_workers, number_updaters, data_generator, worker_function, update_function):
    work_processes = {}
    update_processes = {}
    lock = mp.Lock()
    cores = mp.cpu_count()
    q = mp.Queue(maxsize=cores * 10)
    retq = mp.Queue(maxsize=cores * 10)
    matq = mp.Queue(maxsize=cores * 10)

    print('start workers')
    for i in range(number_workers):
        t = mp.Process(target=worker_function, args=(i, q, retq, matq, lock))
        t.daemon = True
        t.start()
        work_processes[i] = t

    print('start savers')

    for i in range(number_updaters):
        t = mp.Process(target=update_function, args=(i, q, retq, matq, lock))
        t.daemon = True
        t.start()
        update_processes[i] = t

    count = 0
    start = time.time()

    while True:
        time.sleep(.1)
        try:
            data = next(data_generator)
            q.put(data)
            count += 1

        except StopIteration:
            print('stop iteration')

            for p in range(number_workers):
                q.put(None)
            for p in range(number_updaters):
                retq.put(None)

            for p in work_processes:
                work_processes[p].join()
            for p in update_processes:
                update_processes[p].join()

            for p in work_processes:
                work_processes[p].terminate()
            for p in update_processes:
                update_processes[p].terminate()

            break

    gc.collect()
    print('DONE!!!!!')