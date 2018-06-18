from tables import *
import functools
import gc
import multiprocessing as mp
import pandas as pd
import time
import h5sparse
import pickle
from datasketch import MinHashLSH , MinHashLSHForest
from scipy import sparse
from datetime import datetime
import h5py

from pyoma.browser import db

from utils import files_utils, config_utils, pyhamutils, hashutils


class LSHBuilder:

    def __init__(self, h5_oma, saving_path, hog_level=None , numperm = 128):
        print("starting LSH BUILDER with hog level {}".format(hog_level))
        self.h5OMA = h5_oma
        self.db_obj = db.Database(h5_oma)
        self.oma_id_obj = db.OmaIdMapper(self.db_obj)
        self.tree = files_utils.get_tree(self.h5OMA)
        self.taxaIndex, self.reverse = files_utils.generate_taxa_index(self.h5OMA)
        self.numperm = numperm
        self.saving_path = saving_path
        self.datetime = datetime
        self.date_string = "{:%B_%d_%Y_%H_%M}".format(datetime.now())

        # define functions
        self.HAM_PIPELINE = functools.partial(pyhamutils.get_ham_treemap_from_row, tree=self.tree)
        self.HASH_PIPELINE = functools.partial(hashutils.tree2hashes_from_row, events=['duplication', 'gain', 'loss', 'presence'], combination=True , nperm =numperm )
        self.ROW_PIPELINE = functools.partial(hashutils.tree2mat, taxaIndex=self.taxaIndex)
        self.READ_ORTHO = functools.partial(pyhamutils.get_orthoxml, db_obj=self.db_obj)
        if hog_level is not None:
            self.allowed_families = files_utils.get_allowed_families(self.db_obj, hog_level)
            print(self.allowed_families)
            print(len(self.allowed_families))
            self.columns = len(self.allowed_families)
            self.rows = len(self.allowed_families)
        else:
            self.allowed_families = None
            self.columns = len(self.taxaIndex)
            self.rows = len(self.h5OMA.root.OrthoXML.Index)

    def generates_dataframes(self, size=100):

        families = {}
        for i, row in enumerate(self.h5OMA.root.OrthoXML.Index):
            fam = row[0]
            if self.allowed_families is None or fam in self.allowed_families:
                families[fam] = {'ortho': self.READ_ORTHO(fam)}
                if len(families) > size:
                    pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                    pd_dataframe['Fam'] = pd_dataframe.index
                    families = {}
                    yield pd_dataframe

    def worker(self, i, q, retq, matq, l):

        print('worker init ' + str(i))
        while True:
            df = q.get()
            if df is not None:

                df['tree'] = df[['Fam', 'ortho']].apply(self.HAM_PIPELINE, axis=1)
                df['hash'] = df[['Fam', 'tree']].apply(self.HASH_PIPELINE, axis=1)

                # TODO uncomment for rows
                # df['rows'] = df['tree'].apply(self.ROW_PIPELINE)
                retq.put(df[['Fam', 'hash']])
                # TODO uncomment for rows
                # matq.put(df['rows'])
            else:
                print('Worker done' + str(i))
                break

    def saver(self, i, q, retq, matq, l):

        print_start = time.clock()
        save_start = time.clock()
        global_time = time.clock()

        chunk_size = 100
        count = 0

        dataset_names = ['duplication', 'gain', 'loss', 'presence']
        threshold = 0.7

        lsh = MinHashLSH(threshold=threshold, num_perm=self.numperm)
        forest = MinHashLSHForest(num_perm=self.numperm)

        with open(self.saving_path + self.date_string + 'errors.txt', 'a') as hashes_error_files:
            with h5py.File(self.saving_path + self.date_string + 'hashes.h5', 'w', libver='latest') as h5hashes:
                datasets = {}
                for dataset_name in dataset_names:
                    print(dataset_name)
                    if dataset_name not in list(h5hashes.keys()):
                        h5hashes.create_dataset(dataset_name, (chunk_size, 0), maxshape=(None, None), dtype='int32')
                    datasets[dataset_name] = h5hashes[dataset_name]
                print(datasets)
                print('saver init ' + str(i))

                while True:
                    this_dataframe = retq.get()

                    print(time.clock() - global_time)
                    print(count)

                    if this_dataframe is not None:

                        hashes = this_dataframe['hash'].to_dict()

                        for fam in hashes:

                            if hashes[fam] is not None:
                                for famhashname in hashes[fam]['dict']:
                                    lsh.insert(famhashname, hashes[fam]['dict'][famhashname])
                                    forest.add(famhashname, hashes[fam]['dict'][famhashname])
                                hashvals = hashes[fam]['hashes']
                                for event in hashvals:
                                    if len(datasets[event]) < fam + 10:
                                        datasets[event].resize((fam + chunk_size, len(hashvals[event].hashvalues)))
                                    datasets[event][fam, :] = hashvals[event].hashvalues
                                h5hashes.flush()

                            else:
                                print('error')
                                print(fam)
                                hashes_error_files.write(str(fam) + '\n')
                        if time.clock() - print_start > 60:
                            print(this_dataframe['Fam'].max())
                            print(time.clock() - global_time)
                            print_start = time.clock()

                        if time.clock() - save_start > 2000:
                            print(this_dataframe['Fam'].max())
                            print('saving')
                            with open(self.saving_path + self.date_string + '_' + str(threshold) + '_' + 'newlsh.pkl', 'wb') as lsh_out:
                                pickle.dump(lsh, lsh_out, -1)
                            with open(self.saving_path + self.date_string + 'newlshforest.pkl', 'wb') as forestout:
                                pickle.dump(forest, forestout, -1)
                            save_start = time.clock()
                            print(time.clock() - global_time)
                        count += len(this_dataframe)
                    else:
                        with open(self.saving_path + self.date_string + '_' + str(threshold) + '_' + 'newlsh.pkl',
                                  'wb') as lsh_out:
                            pickle.dump(lsh, lsh_out, -1)
                        with open(self.saving_path + self.date_string + 'newlshforest.pkl', 'wb') as forestout:
                            pickle.dump(forest, forestout, -1)
                        print('DONE UPDATER' + str(i))
                        break

    def run_pipeline(self):
        self.mp_with_timeout(number_workers=int(mp.cpu_count() / 1.5), number_updaters=1,
                             data_generator=self.generates_dataframes(100), worker_function=self.worker,
                             update_function=self.saver)

    def matrix_updater(self, i, q, retq, matq, l, rows, columns):
        hog_mat = sparse.csr_matrix((rows, columns))
        save_start = time.clock()

        while True:
            rows = matq.get()
            if rows is not None:
                for fam, sparse_row in rows.itesm():
                    hog_mat[fam, :] = sparse_row
                if time.clock() - save_start > 2000:
                    with h5sparse.File(self.saving_path + self.date_string + "matrix.h5", 'w') as h5matrix:
                        h5matrix.create_dataset('hogmat', data=hog_mat)

            else:

                with h5sparse.File(self.saving_path + self.date_string + "matrix.h5", 'w') as h5matrix:
                    h5matrix.create_dataset('hogmat', data=hog_mat)
                print('DONE MAT UPDATER' + str(i))

                break

    @staticmethod
    def mp_with_timeout(number_workers, number_updaters, data_generator, worker_function, update_function ):
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

        while True:
            # time.sleep(.1)
            try:
                data = next(data_generator)
                q.put(data)
                count += 1

            except StopIteration:
                print('stop iteration')

                for p in range(2*number_workers):
                    q.put(None)

                # for p in work_processes:
                #     work_processes[p].join()

                break

        for p in range(2*number_updaters):
            retq.put(None)

        for p in update_processes:
            update_processes[p].join()

        # for p in work_processes:
        #     work_processes[p].terminate()
        # for p in update_processes:
        #     update_processes[p].terminate()

        gc.collect()
        print('DONE!')


if __name__ == '__main__':

    # hyper params
    num_perm = config_utils.num_perm

    # for use with weighted minhash functions
    # overall weight of category
    lossweight = [0,1]
    presencweight = [0,1]
    gainweight = [0,1]
    dupweight=[0,1]

    # bleed to neighbors up and down
    lossbleed = [0,1]
    presencebleed = [0,1]
    gainbleed = [0,1]
    dupbleed=[0,1]

    # importance given to taxonomic levels
    losslin = [0,1]
    presencelin = [0,1]
    gainlin = [0,1]
    duplin=[0,1]

    with open_file(config_utils.omadir + 'OmaServer.h5', mode="r") as h5OMA:
        # loop with bayes opt over hyper params

        # build lsh
        lsh_builder = LSHBuilder(h5_oma=h5OMA, saving_path=config_utils.datadir, numperm=num_perm)
        lsh_builder.run_pipeline()

        # run validation

        # output score and params
