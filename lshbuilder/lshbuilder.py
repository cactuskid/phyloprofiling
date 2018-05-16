from tables import *
import functools
import gc
import multiprocessing as mp
import pandas as pd
import time
import h5sparse
import pickle
from datasketch import MinHashLSH, MinHashLSHForest
from scipy import sparse
from datetime import datetime
import h5py

from pyoma.browser import db

from utils import pyhamutils
from utils import hashutils
from utils import files_utils
from utils import config_utils


class LSHBuilder:

    def __init__(self, h5_oma, saving_path):
        self.h5OMA = h5_oma
        self.db_obj = db.Database(h5_oma)
        self.oma_id_obj = db.OmaIdMapper(self.db_obj)
        self.dic, self.tree = files_utils.get_species_tree_replacement_dic(h5_oma, self.oma_id_obj)
        self.taxaIndex, self.reverse = files_utils.generate_taxa_index(self.tree)

        self.saving_path = saving_path
        self.datetime = datetime
        self.date_string = "{:%B_%d_%Y_%H_%M}".format(datetime.now())

        # define functions
        self.HAM_PIPELINE = functools.partial(pyhamutils.get_ham_treemap, species_tree=self.tree)
        self.HASH_PIPELINE = functools.partial(hashutils.tree2hashes_from_row,
                                               events=['duplication', 'gain', 'loss', 'presence'], combination=True)
        self.ROW_PIPELINE = functools.partial(hashutils.tree2mat, taxaIndex=self.taxaIndex)
        self.READ_ORTHO = functools.partial(pyhamutils.get_orthoxml, dbObj=self.db_obj, species_tree=self.tree,
                                            replacement_dic=self.dic)

        self.columns = len(self.taxaIndex)
        self.rows = len(self.h5OMA.root.OrthoXML.Index)

    def generates_dataframes(self, size=100):

        families = {}
        for i, row in enumerate(self.h5OMA.root.OrthoXML.Index):
            fam = row[0]
            families[fam] = {'ortho': self.READ_ORTHO(fam)}
            if len(families) > size:
                pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                pd_dataframe['Fam'] = pd_dataframe.index
                families = {}
                yield pd_dataframe

    def worker(self, i, q, retq, matq):

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

        dataset_names = ['Fam', 'duplication', 'gain', 'loss', 'presence']
        threshold = 0.7

        lsh = MinHashLSH(threshold=threshold, num_perm=128)
        forest = MinHashLSHForest(num_perm=128)

        with open(self.saving_path + self.date_string + 'errors.txt', 'w') as hashes_error_files:
            with h5py.File(self.saving_path + self.date_string + 'hashes.h5', 'w', libver='latest') as h5hashes:
                datasets = {}
                for dataset_name in dataset_names:
                    print(dataset_name)
                    if dataset_name not in list(h5hashes.keys()):
                        # TODO why this one is useless ??
                        # dataset = h5hashes.create_dataset(dataset_name, (chunk_size, 0), maxshape=(None, None),
                        #                                  dtype='int32')
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

                                for i, event in enumerate(hashvals):
                                    dataset = dataset_names[i]
                                    if len(datasets[dataset]) < fam + 10:
                                        datasets[dataset].resize((fam + chunk_size, len(hashvals[event])))
                                    datasets[dataset][fam, :] = hashvals[event]
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
                            with open(self.saving_path + self.date_string + '_' + str(threshold) + '_' + 'newlsh.pkl',
                                      'wb') as lsh_out:
                                with open(self.saving_path + self.date_string + 'newlshforest.pkl', 'wb') as forestout:
                                    pickle.dump(lsh, lsh_out, -1)
                                    pickle.dump(forest, forestout, -1)
                            save_start = time.clock()
                            print(time.clock() - global_time)
                        count += len(this_dataframe)
                    else:
                        with open(self.saving_path + self.date_string + '_' + str(threshold) + '_' + 'newlsh.pkl',
                                  'wb') as lsh_out:
                            with open(self.saving_path + self.date_string + 'newlshforest.pkl', 'wb') as forestout:
                                pickle.dump(lsh, lsh_out, -1)
                                pickle.dump(forest.index(), forestout, -1)

                        print('DONE UPDATER' + str(i))
                        break

    def run_pipeline(self):
        self.mp_with_timeout(number_workers=int(mp.cpu_count() / 2), number_updaters=1,
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


if __name__ == '__main__':

    with open_file(config_utils.omadir + 'OmaServer.h5', mode="r") as h5OMA:

        lsh_builder = LSHBuilder(h5_oma=h5OMA, saving_path="saving path")
        lsh_builder.run_pipeline()
