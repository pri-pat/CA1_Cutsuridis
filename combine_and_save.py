#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 14:23:25 2021

@author: daisyl
"""
import numpy as np

def combine_and_save(fname):
    try:
        numpatt = int(str(fname)[-2:])
    except:
        try: 
            numpatt = int(str(fname)[-1])
        except:
            numpatt = int(str(fname[-1 - len('combined')]))
    #open and read memory
    realmem = np.loadtxt("Weights/patts" + fname + ".dat")
   
    comb_mem = realmem.copy()
      
    chunklen = 100//numpatt
        
    for i in range(numpatt):
        if i == 0:
            comb_mem[:,i] = realmem[:,i].copy() 
        else:
            start = i* chunklen
            stop = (i + 1)* chunklen
            (comb_mem[:, 0])[start:stop] = (realmem[:,i])[start:stop]
            for n in range(100):
                comb_mem[n, i] = 0
 
    return comb_mem