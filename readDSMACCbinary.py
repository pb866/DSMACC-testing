import os, multiprocessing
import pandas as pd
import numpy as np 
from scipy.io import FortranFile
import matplotlib.pyplot as plt



def read_fbin(filename):
    ''' this reads each written binary instance itteratively'''
    f = FortranFile(filename, 'r')
    array = []
    while True:
        try: 
            array.append(f.read_reals(dtype=np.float_))
        except TypeError: 
            break
    #array = np.reshape(array, (nspecs,-1))
    
    f.close()
    return array
    
 
global nspecs    

species = ''.join(tuple(open('spec.names'))).replace('\n','').replace(' ','')[:-1].split(',') 
nspecs = len(species)-1
array = read_fbin('1_def.spec')


specs = pd.DataFrame(array , columns = species)

