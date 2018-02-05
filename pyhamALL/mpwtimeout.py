
import pickle 
import multiprocessing as mp
import time
import gc
from multiprocessing import TimeoutError, Manager

import signal
from contextlib import contextmanager


from dask import dataframe
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


    

def updatefunction(startobject, update_data):
    #for a dictionary
 
    try:
        startobject.update(update_data)
    except:
        print( update_data)
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



def daskupdater(i,q,retq,l , updatefunction ,startobject, saveobject ):
    #build up a dictionary
    #once enoguh entries are there transform it into a pandas dataframe
    #shove it into the final Dask dataframe and update it
    pnumber = i
    while True:
        time.sleep(.1)
        update_data = retq.get()
        if update_data == 'DONE':
            break
        elif update_data == 'SAVE':
            
            df = pd.from_dict()


            startobject = {}
        else:
            startobject = updatefunction(startobject , update_data)

            
def mp_with_timeout(nworkers, nupdaters, startobject , saveobject , datagenerator , workerfunction, updatefunction, timeout = 60, saveinterval = 600  ):

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
        t = mp.Process(target=updater, args=(i,q,retq,l , updatefunction ,startobject , saveobject+str(i)+'.pkl'  ) ) 
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
        if count % 100 == 0 :
            retq.put('SAVE')
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
