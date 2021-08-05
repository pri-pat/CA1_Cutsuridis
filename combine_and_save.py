#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 14:23:25 2021

@author: daisyl
"""
import numpy as np

def combine_and_save(fname):
    numpatt = int(str(fname)[-1])
    realmem = ""
    
    #open and read memory
    realmem = np.loadtxt("Weights/patts" + fname + ".dat")
   
    comb_mem = realmem.copy()
      
    chunklen = 100//numpatt
    
    for i in range(numpatt):
        if i == 0:
            comb_mem[:,0] = realmem[:,0].copy() 
        else:
            start = i* chunklen
            stop = (i + 1)* chunklen
            (comb_mem[:, 0])[start:stop] = (realmem[:,i])[start:stop]
    
    return comb_mem