from tables import *
import functools
import gc
import multiprocessing as mp
import pandas as pd
import time
import pickle
from datasketch import MinHashLSH, MinHashLSHForest
from scipy import sparse
from datetime import datetime
import h5py
import redis

from pyoma.browser import db

from utils import files_utils, config_utils, pyhamutils, hashutils


class LSHBuilder:

    def __init__(self, h5_oma, saving_folder , saving_name=None , numperm = 128, tax_filter=None):
        self.h5OMA = h5_oma
        self.db_obj = db.Database(h5_oma)
        self.oma_id_obj = db.OmaIdMapper(self.db_obj)
        self.tax_filter = tax_filter
        self.tree_string, self.tree_ete3 = files_utils.get_tree(self.h5OMA)
        self.taxaIndex, self.reverse = files_utils.generate_taxa_index(self.tree_ete3)

        self.numperm = numperm
        self.saving_folder = saving_folder
        self.datetime = datetime
        self.date_string = "{:%B_%d_%Y_%H_%M}".format(datetime.now())
        if saving_name:
            self.saving_path =self.saving_folder + saving_name
        else:
            self.saving_path = self.saving_folder + self.date_string

        # define functions
        self.HAM_PIPELINE = functools.partial(pyhamutils.get_ham_treemap_from_row, tree=self.tree_string , tax_filter = self.tax_filter)
        self.HASH_PIPELINE = functools.partial(hashutils.tree2hashes_from_row, events=['duplication', 'gain', 'loss', 'presence'], combination=True , nperm =numperm )
        self.ROW_PIPELINE = functools.partial(hashutils.tree2mat, taxa_index=self.taxaIndex)
        self.READ_ORTHO = functools.partial(pyhamutils.get_orthoxml, db_obj=self.db_obj)

        self.columns = len(self.taxaIndex)
        self.rows = len(self.h5OMA.root.OrthoXML.Index)

    def generates_dataframes(self, size=100, minhog_size=5, maxhog_size=None):
        families = {}
        for i, row in enumerate(self.h5OMA.root.OrthoXML.Index):
            fam = row[0]
            ## TODO: add further quality check here for hog_size / hogspread
            ortho_fam = self.READ_ORTHO(fam)
            hog_size = ortho_fam.count('<species name=')

            if (maxhog_size is None or hog_size < maxhog_size) and (minhog_size is None or hog_size > minhog_size):
                families[fam] = {'ortho': ortho_fam}

            if len(families) > size:
                pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
                pd_dataframe['Fam'] = pd_dataframe.index
                yield pd_dataframe
                families = {}

    def worker(self, i, q, retq, matq, l):

        print('worker init ' + str(i))
        while True:
            df = q.get()
            if df is not None:
                df['tree'] = df[['Fam', 'ortho']].apply(self.HAM_PIPELINE, axis=1)
                df['hash'] = df[['Fam', 'tree']].apply(self.HASH_PIPELINE, axis=1)
                #df['rows'] = df['tree'].apply(self.ROW_PIPELINE)
                retq.put(df[['Fam', 'hash']])
                #matq.put(df[['Fam', 'rows']])
            else:
                print('Worker done' + str(i))
                break

    def saver(self, i, q, retq, matq, l):
        print_start = time.clock()
        save_start = time.clock()
        global_time = time.clock()
        chunk_size = 100
        count = 0

        dataset_names = ['duplication', 'gain', 'loss', 'presence', 'weighted']
        threshold = 0.7
        if config_utils.clear_redisLSH == True:
            #flush the LSH DB
            r = redis.StrictRedis(host='10.0.63.33', port=6379, db=2)
            r.flushdb()

        lsh = MinHashLSH(
            threshold=threshold,
            num_perm=self.numperm,
            storage_config={'type': 'redis', 'redis': {'host': '10.0.63.33', 'port': 6379, 'db': 2}})

        forest = MinHashLSHForest(num_perm=self.numperm)

        #create datasets
        if self.tax_filter is not None:
            taxstr = 'NoFilter'
        else:
            taxstr = str(self.tax_filter)

        with open(self.saving_path + 'errors.txt', 'w') as hashes_error_files:
            with h5py.File(self.saving_path + 'hashes.h5', 'a', libver='latest') as h5hashes:
                datasets = {}
                for dataset_name in dataset_names:
                    if dataset_name+'_'+taxstr not in list(h5hashes.keys()):
                        print('creating dataset')
                        print(dataset_name)
                        print('filtered at taxonomic level:'+taxstr)
                        h5hashes.create_dataset(dataset_name+'_'+taxstr, (chunk_size, 0), maxshape=(None, None), dtype='int32')
                    datasets[dataset_name] = h5hashes[dataset_name+'_'+taxstr]

                print(datasets)
                print('saver init ' + str(i))
                while True:
                    this_dataframe = retq.get()
                    print(str(time.clock() - global_time)+'seconds elapsed')
                    print(str(this_dataframe.Fam.max())+ 'fam num')
                    print(str(count) + 'hogs done')
                    if this_dataframe is not None:
                        hashes = this_dataframe['hash'].to_dict()
                        for fam in hashes:
                            if hashes[fam] is not None:
                                for famhashname in hashes[fam]['dict']:
                                    #lsh.insert(famhashname, hashes[fam]['dict'][famhashname])
                                    forest.add(famhashname, hashes[fam]['dict'][famhashname])
                                hashvals = hashes[fam]['hashes']
                                for event in hashvals:
                                    if len(datasets[event]) < fam + 10:
                                        datasets[event].resize((fam + chunk_size, len(hashvals[event].hashvalues)))
                                    #set the hasvalues for each event and for the weighted minhash
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

                        if time.clock() - save_start > 1000:
                            print(this_dataframe['Fam'].max())
                            print('saving')
                        #     with open(self.saving_path + self.date_string + '_' + str(threshold) + '_' + 'newlsh.pkl', 'wb') as lsh_out:
                        #         pickle.dump(lsh, lsh_out, -1)
                            with open(self.saving_path + 'newlshforest.pkl', 'wb') as forestout:
                                pickle.dump(forest, forestout, -1)
                            save_start = time.clock()
                            print(time.clock() - global_time)
                        count += len(this_dataframe)
                    else:
                        # with open(self.saving_path + self.date_string + '_' + str(threshold) + '_' + 'newlsh.pkl',
                        #           'wb') as lsh_out:
                        #     pickle.dump(lsh, lsh_out, -1)
                        with open(self.saving_path + 'newlshforest.pkl', 'wb') as forestout:
                            pickle.dump(forest, forestout, -1)
                        print('DONE UPDATER' + str(i))
                        break

    def run_pipeline(self):
        # self.mp_with_timeout(number_workers=int(mp.cpu_count() / 1.5), number_updaters=1,
        #                      data_generator=self.generates_dataframes(100), worker_function=self.worker,
        #                      update_function=self.saver)

        functype_dict = {'worker': (self.worker, int(mp.cpu_count()/2), True), 'updater': (self.saver, 1, False),
                         'matrix_updater': (self.matrix_updater, 1, False)}

        self.mp_with_timeout(functypes=functype_dict, data_generator=self.generates_dataframes(100))

    def matrix_updater(self, i, q, retq, matq, l):
        hog_mat = None
        save_start = time.clock()
        print('hogmat saver init ' + str(i))
        while True:
            rows = matq.get()
            if rows != None:
                for index, row in rows.iterrows():
                    if row is not None:
                        sparse_row = row['rows']
                        fam = int(row['Fam'])
                        try:
                            if not hog_mat:
                                hog_mat = sparse.csr_matrix((fam+10000, sparse_row.shape[1]))
                        except ValueError:
                            pass
                        if hog_mat.shape[0] < fam:
                            print('extend HOGMAT')
                            num_rows_to_add = fam - hog_mat.shape[0] + 10000
                            new_hog_mat = sparse.csr_matrix((num_rows_to_add, sparse_row.shape[1]))
                            hog_mat = sparse.vstack([hog_mat, new_hog_mat])
                        hog_mat[fam, :] = sparse_row
                    else:
                        break

                    if time.clock() - save_start > 2000:
                        # with h5sparse.File(self.saving_path + self.date_string + "matrix.h5", 'w') as h5matrix:
                            # h5matrix.create_dataset('hogmat', data=hog_mat)
                        print('saving HOGMAT')
                        with open(self.saving_path + self.date_string +'_matnum_'+ str(i) + "matrix.pkl", 'wb') as handle:
                            pickle.dump(hog_mat, handle, -1)
                        save_start = time.clock()
            else:
                break

        with open(self.saving_path + '_matnum_'+ str(i) + "matrix.pkl", 'wb') as handle:
            pickle.dump(hog_mat, handle, -1)
        print('DONE MAT UPDATER' + str(i))


    @staticmethod
    def mp_with_timeout(functypes, data_generator):
        work_processes = {}
        update_processes = {}
        lock = mp.Lock()
        cores = mp.cpu_count()
        q = mp.Queue(maxsize=cores * 10)
        retq = mp.Queue(maxsize=cores * 10)
        matq = mp.Queue(maxsize=cores * 10)

        work_processes = {}
        print('start workers')
        for key in functypes:
            worker_function, number_workers, joinval = functypes[key]
            work_processes[key] = []
            for i in range(int(number_workers)):
                t = mp.Process(target=worker_function, args=(i, q, retq, matq, lock))
                t.daemon = True
                work_processes[key].append(t)
        for key in work_processes:
            for process in work_processes[key]:
                process.start()

        count = 0

        while True:
            # time.sleep(.1)
            try:
                data = next(data_generator)
                q.put(data)
                count += 1

            except StopIteration:
                print('stop iteration')
                for key in work_processes:
                    for _ in work_processes[key]:
                        q.put(None)
                break

        for key in work_processes:
            worker_function, number_workers , joinval = functypes[key]
            if joinval == True:
                for process in work_processes[key]:
                    process.join()

        for key in work_processes:
            worker_function, number_workers, joinval = functypes[key]
            if joinval == False:
                for _ in work_processes[key]:
                    retq.put(None)
                    matq.put(None)

        gc.collect()
        print('DONE!')


if __name__ == '__main__':

    # hyper params
    num_perm = config_utils.num_perm

    # for use with weighted minhash functions
    # overall weight of category
    lossweight = [0, 1]
    presencweight = [0, 1]
    gainweight = [0, 1]
    dupweight=[0, 1]

    # bleed to neighbors up and down
    lossbleed = [0, 1]
    presencebleed = [0, 1]
    gainbleed = [0, 1]
    dupbleed=[0, 1]

    # importance given to taxonomic levels
    losslin = [0, 1]
    presencelin = [0, 1]
    gainlin = [0, 1]
    duplin=[0, 1]

    with open_file(config_utils.omadir + 'OmaServer.h5', mode="r") as h5OMA:
        # loop with bayes opt over hyper params

        # build lsh
        lsh_builder = LSHBuilder(h5_oma=h5OMA, saving_folder =config_utils.datadir , saving_name = 'FinalDB', numperm=256, tax_filter=2759)
        lsh_builder.run_pipeline()

        # run validation

        # output score and params
