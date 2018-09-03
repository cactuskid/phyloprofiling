from tables import *
import functools
import gc
import multiprocessing as mp
import pandas as pd
import time as t
import pickle
from datasketch import MinHashLSH,   MinHashLSHForest
from scipy import sparse
from datetime import datetime
import h5py
import redis

from pyoma.browser import db

from utils import files_utils, config_utils, pyhamutils, hashutils
from datasketch import WeightedMinHashGenerator


class LSHBuilder:

    def __init__(self, h5_oma, saving_folder , saving_name=None , numperm = 128,  treeweights= None , taxfilter = None, taxmask= None , lambdadict= None, start= None):
        self.h5OMA = h5_oma
        self.db_obj = db.Database(h5_oma)
        self.oma_id_obj = db.OmaIdMapper(self.db_obj)
        self.tax_filter = taxfilter
        self.tax_mask = taxmask
        self.tree_string, self.tree_ete3 = files_utils.get_tree(self.h5OMA)
        self.taxaIndex, self.reverse = files_utils.generate_taxa_index(self.tree_ete3)
        self.numperm = numperm

        if treeweights is None:
            self.treeweights = hashutils.generate_treeweights(self.tree_ete3  , self.taxaIndex , taxfilter, taxmask , lambdadict, start)
        else:
            self.treeweights = treeweights


        self.saving_folder = saving_folder
        self.datetime = datetime
        self.date_string = "{:%B_%d_%Y_%H_%M}".format(datetime.now())
        self.saving_name= saving_name
        if saving_name:
            self.saving_path =self.saving_folder + saving_name
        else:
            self.saving_path = self.saving_folder + self.date_string

        # define minhash generator. taxa stay fixedself.

        wmg = WeightedMinHashGenerator(3*len(self.taxaIndex), sample_size=numperm, seed=1)
        self.wmg = wmg

        self.HAM_PIPELINE = functools.partial(pyhamutils.get_ham_treemap_from_row, tree=self.tree_string , tax_filter = self.tax_filter)

        self.HASH_PIPELINE = functools.partial(hashutils.row2hash , taxaIndex=self.taxaIndex  , treeweights=self.treeweights , wmg=wmg )

        self.READ_ORTHO = functools.partial(pyhamutils.get_orthoxml, db_obj=self.db_obj)

        self.columns = len(self.taxaIndex)
        self.rows = len(self.h5OMA.root.OrthoXML.Index)

    def generates_dataframes(self, size=100, minhog_size=None, maxhog_size=None):
        families = {}
        start = -1

        for i, row in enumerate(self.h5OMA.root.OrthoXML.Index):
            if i > start:
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

        pd_dataframe = pd.DataFrame.from_dict(families, orient='index')
        pd_dataframe['Fam'] = pd_dataframe.index
        yield pd_dataframe
        print('last dataframe sent')
        families = {}


    def worker(self, i, q, retq, matq, l):

        print('worker init ' + str(i))
        while True:
            df = q.get()
            if df is not None :
                df['tree'] = df[['Fam', 'ortho']].apply(self.HAM_PIPELINE, axis=1)
                df[['hash','rows']] = df[['Fam', 'tree']].apply(self.HASH_PIPELINE, axis=1)
                retq.put(df[['Fam', 'hash']])
                #matq.put(df[['Fam', 'rows']])
            else:
                print('Worker done' + str(i))
                break


    def saver(self, i, q, retq, matq, l):
        print_start = t.time()
        save_start = t.time()
        global_time = t.time()
        chunk_size = 100
        count = 0
        threshold = 0.98

        if config_utils.clear_redisLSH == True:
            #flush the LSH DB
            r = redis.StrictRedis(host='10.0.63.33', port=6379, db=2)
            r.flushdb()

        lsh = MinHashLSH(
            threshold=threshold,
            num_perm=self.numperm
            )#storage_config={'type': 'redis', 'redis': {'host': '10.0.63.33', 'port': 6379, 'db': 2}})

        forest = MinHashLSHForest(num_perm=self.numperm)

        #create datasets
        if self.tax_filter is None:
            taxstr = 'NoFilter'
        if self.tax_mask is None:
            taxstr+= 'NoMask'
        else:
            taxstr = str(self.tax_filter)
        dataset_name = self.saving_name+'_'+taxstr
        with open(self.saving_path + 'errors.txt', 'w') as hashes_error_files:
            with h5py.File(self.saving_path + 'hashes.h5', 'w', libver='latest') as h5hashes:
                datasets = {}
                if dataset_name not in h5hashes.keys():
                    print('creating dataset')
                    print(dataset_name)
                    print('filtered at taxonomic level:'+taxstr)
                    h5hashes.create_dataset(dataset_name+'_'+taxstr, (chunk_size, 0), maxshape=(None, None), dtype='int32')
                    datasets[dataset_name] = h5hashes[dataset_name+'_'+taxstr]
                    print(datasets)
                print('saver init ' + str(i))
                while True:
                    this_dataframe = retq.get()

                    if not this_dataframe.empty:
                        hashes = this_dataframe['hash'].to_dict()
                        print(str(t.time() - global_time)+'seconds elapsed')
                        print(str(this_dataframe.Fam.max())+ 'fam num')
                        print(str(count) + ' done')

                        for fam in hashes:
                            if hashes[fam] is not None:

                                lsh.insert(str(fam), hashes[fam])
                                forest.add(str(fam), hashes[fam])
                                if len(datasets[dataset_name]) < fam + 10:
                                    datasets[dataset_name].resize((fam + chunk_size, len(hashes[fam].hashvalues.ravel())))
                                datasets[dataset_name][fam, :] = hashes[fam].hashvalues.ravel()
                                h5hashes.flush()
                            else:
                                print('error')
                                print(fam)
                                hashes_error_files.write(str(fam) + '\n')

                        if t.time() - save_start > 200:

                            print('saving')
                            forest.index()
                            print(forest.query(hashes[fam],100))
                            #print(lsh.query(hashes[fam]))
                            with open(self.saving_path + 'newlsh.pkl','wb') as lsh_out:
                                lsh_out.write(pickle.dumps(lsh, -1))
                            with open(self.saving_path + 'newlshforest.pkl', 'wb') as forestout:
                                forestout.write(pickle.dumps(forest, -1))

                            save_start = t.time()

                            print('save done at' + str(t.time() - global_time))
                        count += len(this_dataframe)
                    else:
                        print('wrap it up')
                        #wrap it up
                        with open(self.saving_path + 'newlsh.pkl','wb') as lsh_out:
                            lsh_out.write(pickle.dumps(lsh, -1))
                        with open(self.saving_path + 'newlshforest.pkl', 'wb') as forestout:
                            forestout.write(pickle.dumps(forest, -1))

                        print('DONE UPDATER' + str(i))
                        break

    def run_pipeline(self):
        ## TODO: return files saved

        functype_dict = {'worker': (self.worker, int(mp.cpu_count()/2), True), 'updater': (self.saver, 1, False),
                         'matrix_updater': (self.matrix_updater, 0, False)}

        self.mp_with_timeout(functypes=functype_dict, data_generator=self.generates_dataframes(100))

    def matrix_updater(self, i, q, retq, matq, l):
        hog_mat =  sparse.lil_matrix((600000, len(self.taxaIndex)*3))
        save_start = t.time()

        print('hogmat saver init ' + str(i))

        while True:
            rows = matq.get()
            if rows is not None:
                for index, row in rows.iterrows():
                    if row is not None:
                        sparse_row = row['rows']
                        fam = int(row['Fam'])
                        if hog_mat.shape[0] < fam:
                            print('extend HOGMAT')
                            num_rows_to_add = fam - hog_mat.shape[0] + 1000
                            new_hog_mat = sparse.lil_matrix((num_rows_to_add, len(self.taxaIndex)*3 ))
                            hog_mat = sparse.vstack([hog_mat, new_hog_mat])
                        hog_mat[fam, :] = sparse_row


                    if t.time() - save_start > 500:
                        # with h5sparse.File(self.saving_path + self.date_string + "matrix.h5", 'w') as h5matrix:
                            # h5matrix.create_dataset('hogmat', data=hog_mat)
                        print('saving HOGMAT')
                        with open(self.saving_path + '_matnum_'+ str(i) + "matrix.pkl", 'wb') as handle:
                            pickle.dump(hog_mat, handle, -1)
                        save_start = t.time()
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
        q = mp.Queue(maxsize=cores * 40)
        retq = mp.Queue(maxsize=cores * 40)
        matq = mp.Queue(maxsize=cores * 40)

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

        for data in data_generator:
            q.put(data)
        print('done spooling data')

        for key in work_processes:
            for i in range(2):
                for _ in work_processes[key]:
                    q.put(None)

        print('joining processes')
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

        for key in work_processes:
            worker_function, number_workers , joinval = functypes[key]
            if joinval == False:
                for process in work_processes[key]:
                    process.join()

        gc.collect()
        print('DONE!')


if __name__ == '__main__':

    # hyper params
    num_perm = config_utils.num_perm
    startdict={'presence':1, 'loss':1, 'dup':1}
    lambdadict={'presence':1, 'loss':1, 'dup':1}

    with open_file(config_utils.omadir + 'OmaServer.h5', mode="r") as h5_oma:

        lsh_builder = LSHBuilder(h5_oma, saving_folder= config_utils.datadir , saving_name='final', numperm = 128,
        treeweights= None , taxfilter = None, taxmask=None , lambdadict= lambdadict, start= startdict)
        lsh_builder.run_pipeline()
