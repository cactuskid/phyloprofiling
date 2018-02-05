
import pickle 
import multiprocessing as mp
import time
import gc
from multiprocessing import TimeoutError, Manager

import signal
from contextlib import contextmanager


from dask import dataframe as ddf
import pandas as pd 



class TimeoutException(Exception): pass

@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


def updatefunction_dict(startobject, update_data):
    #for a dictionary
    try:
        startobject.update(update_data)
    except:
        print(update_data)
    return startobject



def updatefunction_sparsemat(startobject, update_data):
    #for a scipy spares matrix
    #finish me
    try:
        pass
    except:
        print(update_data)
    return startobject

def worker(i,q,retq,l, timecards , workerfunction  ):
    pnumber = i
    timecards[ pnumber ] = time.time()
    while True:
        timecards[ pnumber ] = time.time()
        time.sleep(.1)
        data = q.get()
        completed = False
        if data == 'DONE':
            break 
        try:
            with time_limit(20):
                update_data = workerfunction(data , l)
                completed = True
        except TimeoutException as e:
            print("Timed out!")
        
        gc.collect()
        if completed == True:
            retq.put(update_data)


def updater(i,q,retq,l , updatefunction ,startobject, saveobject ):
    pnumber = i
    while True:
        time.sleep(.1)
        update_data = retq.get()
        if update_data == 'DONE':
            break
        elif update_data == 'SAVE':
            handle1=open(saveobject , 'wb')
            pickle.dump(startobject, handle1 , -1)
            handle1.close()
        else:
            startobject = updatefunction(startobject , update_data)

def daskupdater(i,q,retq,l , updatefunction ,startobject, saveobject  ):
    #build up a dictionary
    #once enoguh entries are there transform it into a pandas dataframe
    #shove it into the final Dask dataframe and update it
    count = 0
    pnumber = i
    while True:
        time.sleep(.1)
        update_data = retq.get()
        if update_data == 'DONE':
            break
        elif update_data == 'SAVE':
            df = pd.from_dict( startobject )
            startobject = {}
            if count == 0:
                DDF = dask.dataframe.from_pandas(df) 
            else:
                DDF.add(df)
            count += 1
        elif update_data == 'DDFSAVE':
            #save dask dataframe as 
            dask.dataframe.to_hdf5(DDF, saveobject)
        else:
            startobject = updatefunction(startobject , update_data)
    print('final save')
    dask.dataframe.to_hdf5(DDF, saveobject)
    print('daskupdater:DONE')


def hdf5updater( i,q,retq,l , updatefunction ,startobject, saveobject, sparse = True ):
    #output vectors to an hdf5 file
    #finish me
    

    f = tb.open_file(saveobject, 'w')

    if sparse == True:
        filters = tb.Filters(complevel=5, complib='blosc')
        out_data = f.create_earray(f.root, 'data', tb.Float32Atom(), shape=(0,), filters=filters)
        out_ri = f.create_earray(f.root, 'ri', tb.Float32Atom(),shape=(0,), filters=filters)
        out_ci = f.create_earray(f.root, 'ci', tb.Float32Atom(), shape=(0,), filters=filters)
        bl = 1000 
    #if sparse store index
    else:
        bl = 1000

    if sparse == True:
    

        while True:
        time.sleep(.1)
        update_data = retq.get()
        if update_data == 'DONE':
            break
        
        elif update_data == 'SAVE':
            df = pd.from_dict( startobject )
            startobject = {}
            if count == 0:
                DDF = dask.dataframe.from_pandas(df) 
            else:
                DDF.add(df)
            count += 1
        elif update_data == 'DDFSAVE':
            #save dask dataframe as 
            dask.dataframe.to_hdf5(DDF, saveobject)
        else:
            startobject = updatefunction(startobject , update_data)
            
        for i in range(0, l, bl):
          res = retmat[:,i:min(i+bl, l)]
          vals = res.data
          ri, ci = res.nonzero()
          out_data.append(vals)
          out_ri.append(ri)
          out_ci.append(ci)
    
    else:


def mp_with_timeout(nworkers, nupdaters, startobject , saveobject , datagenerator , workerfunction, updatefunction, timeout = 60, saveinterval = 300 , bigsaveinterval = 1000 ):

    wprocesses ={}
    uprocesses ={}    
    l = mp.Lock()
    cores = mp.cpu_count()
    q = mp.Queue( maxsize = cores*10 )
    retq = mp.Queue( maxsize = cores*10 )
    
    #30sec timetout for workers
    manager = Manager()
    timecards = manager.dict()

    for i in range(nworkers):
        t = mp.Process(target=worker, args=(i,q,retq,l, timecards , workerfunction )  ) 
        t.daemon = True
        t.start()
        wprocesses[i] = t


    for i in range(nupdaters):
        t = mp.Process(target=daskupdater, args=(i,q,retq,l , updatefunction ,startobject , saveobject+str(i)+'.pkl'  ) ) 
        t.daemon = True
        t.start()
        uprocesses[i]= t

    count =0
    data = next(datagenerator)
    start = time.time()
    
    while True:
        time.sleep(.1)
        #check for non responsive workers
        #save every so often
        if count % saveinterval == 0 :
            retq.put('SAVE')
        if count % bigsaveinterval == 0 :
            retq.put('DDFSAVE')
        try:
            data = next(datagenerator)
            q.put(data)
            count += 1
        except StopIteration:
            
            print( 'stop iteration')
            for p in range(nworkers):
                q.put('DONE')
                
            
            for p in range(nupdaters):
                retq.put('DONE')
            break
        if count % 100 == 0:
            print( count)
            retq.put('SAVE')
    
    for p in wprocesses:
        wprocesses[p].terminate()
    gc.collect()
    for p in uprocesses:
        uprocesses[p].terminate()
    gc.collect()
    print( 'DONE!!!!!')
