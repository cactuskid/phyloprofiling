import functools
import gc
import multiprocessing as mp
import pickle
import time
from datetime import datetime

import config
import format_files
import h5py
import h5sparse
import pandas as pd
import profileGen
import pyhamPipeline
from datasketch import MinHashLSH, MinHashLSHForest
from pyoma.browser import db
from tables import *

if __name__ == '__main__':

    with open_file(config.omadir + 'OmaServer.h5', mode="r") as h5OMA:
        dbObj = db.Database(h5OMA)
        omaIdObj = db.OmaIdMapper(dbObj)
        dic, tree = format_files.create_species_tree(h5OMA, omaIdObj)
        taxaIndex, reverse = profileGen.generateTaxaIndex(tree)

        HAMPIPELINE = functools.partial(pyhamPipeline.runpyham, species_tree=tree)
        HASHPIPEline = functools.partial(profileGen.DFTree2Hashes)
        ROWPIPELINE = functools.partial(profileGen.Tree2mat, taxaIndex=taxaIndex)
        columns = len(taxaIndex)
        rows = len(h5OMA.root.OrthoXML.Index)

    def genDFs( h5OMA, size = 100):
        fams = {}
        with open_file(config.omadir + 'OmaServer.h5', mode="r") as h5OMA:
            dbObj = db.Database(h5OMA)
            omaIdObj = db.OmaIdMapper(dbObj)
            dic, tree = format_files.create_species_tree(h5OMA, omaIdObj)
            taxaIndex,reverse  = profileGen.generateTaxaIndex(tree)
            READORTHO = functools.partial( pyhamPipeline.readortho , dbObj= dbObj , species_tree=tree , replacement_dic= dic )
            for i,row in enumerate( h5OMA.root.OrthoXML.Index):
                fam = row[0]
                fams[fam] = { 'ortho':READORTHO( fam ) }
                if len(fams)>size:
                    pddf = pd.DataFrame.from_dict(fams, orient= 'index' )
                    pddf['Fam'] = pddf.index
                    fams = {}
                    yield pddf

    def worker(i, q, retq, matq, l):
        print('worker init ' + str(i))
        while True:
            df = q.get()
            if df is not None:
                df['tree'] = df[['Fam', 'ortho']].apply(HAMPIPELINE, axis=1)
                df['hash'] = df[['Fam', 'tree']].apply(HASHPIPEline, axis=1)

                #df['rows'] = df['tree'].apply( ROWPIPELINE )
                retq.put(df[['Fam','hash']])
                #matq.put(df['rows'])
            else:
                print('DONE WORKER'+ str(i))
                break

    def matrixupdater(i, q , retq , matq , l , rows , columns):

            hogmat = scipy.sparseCSR( (rows , columns )  )
            printstart = time.clock()
            savestart = time.clock()
            globaltime = time.clock()

            while True:
                rows = matq.get()
                if rows != None:
                    for fam,sparserow in rows.itesm():
                        hogmat[fam,:]= sparserow
                    if time.clock() -savestart > 2000:
                        with h5sparse.File(config.datadir + datestr + "matrix.h5",'w'  ) as h5matrix:
                            h5matrix.create_dataset( 'hogmat', data= hogmat)

                else:

                    with h5sparse.File(config.datadir + datestr + "matrix.h5",'w'  ) as h5matrix:
                            h5matrix.create_dataset( 'hogmat', data= hogmat)
                    print('DONE MAT UPDATER' + str(i))

                    break



    def saver( i,q,retq,matq, l ):

        dataset_names = ['fam', 'duplication', 'gain', 'loss', 'presence']

        printstart = time.clock()
        savestart = time.clock()
        globaltime = time.clock()

        chunksize = 100
        count = 0
        dataset_names = ['Fam', 'duplication', 'gain', 'loss', 'presence']
        threshold=0.7
        lsh = MinHashLSH(threshold=threshold, num_perm=128)
        forest = MinHashLSHForest(num_perm=128)

        dt = datetime
        datestr = "{:%B_%d_%Y_%H_%M}".format(dt.now())

        with open('errors.txt' , 'w')as errfile:

            with h5py.File(config.datadir+ datestr + 'hashes.h5','w', libver = 'latest') as h5hashes:
                dsets = {}
                for dataset_name in dataset_names:
                    print(dataset_name)
                    if dataset_name not in list(h5hashes.keys()):
                        dataset = h5hashes.create_dataset(dataset_name, (chunksize,0), maxshape=(None, None), dtype = 'int32')
                    dsets[dataset_name] = h5hashes[dataset_name]
                print(dsets)

                print('saver init ' + str(i))
                while True:
                    thisdf = retq.get()
                    print(time.clock()- globaltime)
                    print(count)
                    if thisdf is not None:

                        hashes = thisdf['hash'].to_dict()
                        for fam in hashes:
                            if hashes[fam] is not None:

                                for famhashname in hashes[fam]['dict']:
                                    lsh.insert(famhashname , hashes[fam]['dict'][famhashname])
                                    forest.add(famhashname ,hashes[fam]['dict'][famhashname])
                                hashvals = hashes[fam]['hashes']

                                for i,event in enumerate(hashvals):
                                    dataset = dataset_names[i]
                                    if len(dsets[dataset])< fam+10:
                                        dsets[dataset].resize( (fam +chunksize , len(hashvals[event])  ) )
                                    dsets[dataset][fam,:] = hashvals[event]
                                h5hashes.flush()
                            else:
                                print('error')
                                print(fam)
                                errfile.write(str(fam) + '\n')

                        if time.clock() - printstart > 60:
                            print(thisdf['Fam'].max())
                            print( time.clock()- globaltime)
                            printstart = time.clock()

                        if time.clock() -savestart > 2000:
                            print(thisdf['Fam'].max())
                            print('saving')
                            with open(config.datadir + datestr + '_'+str(threshold)+'_' +'newlsh.pkl' , 'wb') as lshout:
                                    with open(config.datadir +datestr+ 'newlshforest.pkl' , 'wb') as forestout:
                                        pickle.dump(lsh, lshout, -1)
                                        pickle.dump(forest, forestout, -1)
                            savestart = time.clock()
                            print( time.clock()- globaltime)
                        count+= len(thisdf)
                    else:
                        with open(config.datadir + datestr + '_'+str(threshold)+'_' +'newlsh.pkl' , 'wb') as lshout:
                            with open(config.datadir +datestr+ 'newlshforest.pkl' , 'wb') as forestout:

                                pickle.dump(lsh, lshout, -1)
                                pickle.dump(forest.index(), forestout, -1)

                        print('DONE UPDATER' + str(i))
                        break

    def mp_with_timeout(nworkers, nupdaters, datagenerator , workerfunction, updatefunction ):
        wprocesses ={}
        uprocesses ={}
        l = mp.Lock()
        cores = mp.cpu_count()
        q = mp.Queue( maxsize = cores*10 )
        retq = mp.Queue( maxsize = cores*10 )
        matq = mp.Queue( maxsize = cores*10 )

        print('start workers')
        for i in range(nworkers):
            t = mp.Process(target=workerfunction, args=(i,q,retq,matq, l)  )
            t.daemon = True
            t.start()
            wprocesses[i] = t

        print('start savers')

        for i in range(nupdaters):
            t = mp.Process(target=updatefunction, args=(i,q,retq,matq ,l) )
            t.daemon = True
            t.start()
            uprocesses[i] = t

        #for i in range(nupdaters):
        #	t = mp.Process(target=matrixupdater, args=(i,q,retq,matq ,l, rows, columns) )
        #	t.daemon = True
        #	t.start()
        #	uprocesses[i]= t

        count =0
        start = time.time()

        while True:
            time.sleep(.1)
            try:
                data = next(datagenerator)
                q.put(data)
                count += 1

            except StopIteration:
                print('stop iteration')

                for p in range(nworkers):
                    q.put(None)
                for p in range(nupdaters):
                    retq.put(None)

                for p in wprocesses:
                    wprocesses[p].join()
                for p in uprocesses:
                    uprocesses[p].join()

                for p in wprocesses:
                    wprocesses[p].terminate()
                for p in uprocesses:
                    uprocesses[p].terminate()

                break

        gc.collect()
        print( 'DONE!!!!!')

    mp_with_timeout(int(mp.cpu_count()/2), 1, genDFs(h5OMA, 100), worker, saver)
