# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 10:30:39 2020

@author: mbezaire
"""


import pickle

runs2include=["par","sim1"]

performance_list=[]
death_list=[]
electrostim_list=[]

for run in runs2include:
    with open('pyresults/'+ run+'.pkl') as f:  # Python 3: open(..., 'rb')
        spikeout, vout, perf, data2save = pickle.load(f)
        performance_list.append(data2save["performance"])
        death_list.append(data2save["percentDeath"])
        electrostim_list.append(data2save["electrostim"])
        
import matplotlib.pyplot as plt

plt.figure()
plt.plot(death_list,performance_list)
plt.xlabel('% Cell Death')
plt.ylabel('Memory Recall Performance (Scale of 0 to 1)')
plt.show()
plt.savefig('Images/memory_recall_v_death.png')