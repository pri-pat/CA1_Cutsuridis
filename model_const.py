# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 13:05:11 2020

@author: mbezaire
"""
from neuron import h

connect_random_low_start_ = 1  # low seed for mcell_ran4_init()


Pcell2Pcell_weight = 0.001
Pcell2Pcell_delay = 1

Bcell2Pcell_weight = 0.02
Bcell2Pcell_delay = 1
Pcell2Bcell_weight = 0.0005 
Pcell2Bcell_delay = 1
Bcell2Bcell_weight = 0.001
Bcell2Bcell_delay = 1
Bcell2BScell_weight = 0.02
Bcell2BScell_delay = 1
Bcell2OLMcell_weight = 0.0
Bcell2OLMcell_delay = 1

AAcell2Pcell_weight = 0.04
AAcell2Pcell_delay = 1
Pcell2AAcell_weight = 0.0005
Pcell2AAcell_delay = 1

BScell2Pcell_weight = 0.002	
BScell2Pcell_delay = 1
BScell2BScell_weight = 0.002	# Parameter not set in original code
BScell2BScell_delay = 1	# Parameter not set in original code
BScell2Pcell_GABAB_weight = 0.0004
BScell2Pcell_delay = 1
Pcell2BScell_weight = 0.0005
Pcell2BScell_delay = 1
BScell2Bcell_weight = 0.01
BScell2Bcell_delay = 1

OLMcell2Pcell_weight = 0.04
OLMcell2Pcell_GABAB_weight = 0.0004
OLMcell2Pcell_delay = 1
OLMcell2Bcell_weight = 0.01
OLMcell2Bcell_delay = 1
Pcell2OLMcell_weight = 0.00005
Pcell2OLMcell_delay = 1
OLMcell2Bcell_weight = 0.0
OLMcell2Bcell_delay = 1

h('STARTDEL = 50')    # msecs
h('THETA = 250')    # msecs (4 Hz)
h('GAMMA = 25')    # msecs (40 Hz)
h('ECCA3DEL = 9')    # msecs
h('SIMDUR = 0')    # msecs

# h.SIMDUR = h.STARTDEL + (h.THETA*2) # Scaled down for testing code for technical bugs only
h.SIMDUR = h.STARTDEL + (h.THETA*h.numCycles)    # simulation duration (msecs)


# Septal inhibition
SEPNUM = 1000    # number of SEP spikes
SEPSTART = h.STARTDEL+(h.THETA/12)    # time of first SEP spike
SEPINT = 20    # SEP spike ISI (during burst)
SEPNOISE = 0.4    # SEP ISI noise
SEPBINT = 2*h.THETA/3    # SEP interburst interval
SEPBLEN = h.THETA/3    # SEP burst length
SEPWGT = 0.02	# SEP weight to BCs and AACs
SEPWGTL = 0.0002	# SEP weight to BSCs and OLMs
SEPDEL = 1	# SEP delay

# Background excitation (not used)
ENUM = 0    # number of spikes
ESTART = 0    # time of first spike
EINT = 200    # spike ISI
ENOISE = 1    # ISI noise
EWGT = 0.001    # excitatory weights (AMPA)
ENWGT = 0.002    # excitatory weights (NMDA)
EDEL = 1    # delay (msecs)

# EC excitation
ECPATT = 1    # index of output pattern
ECNUM = 1000    # number of EC spikes
ECSTART = h.STARTDEL    # time of first EC spike
ECINT = h.GAMMA    # EC spike ISI
ECNOISE = 0.2    # EC ISI noise
ECWGT = 0.0    # EC weight to PCs
#ECWGT = 0.001    # EC weight to PCs
ECDEL = 1    # EC delay
EIWGT = 0.00015    # excitatory weights to INs
EIDEL = 1    # delay (msecs)

# Cue (CA3) excitation
CNUM = 1000    # number of cue spikes
CSTART = h.STARTDEL+h.ECCA3DEL    # time of first cue spike
CINT = h.GAMMA    # cue spike ISI
CNOISE = 0.2    # cue ISI noise
CHWGT = 0.0015    # cue weight
CLWGT = 0.0005    # unlearnt weight (usually 0)
CNWGT = 0.0005    # excitatory weights (NMDA)
CDEL = 1    # cue delay

# STDP configuration
STDPDFAC = 0    # depression factor
STDPPFAC = 0    # potentiation factor
#STDPDFAC = 0.2    # depression factor
#STDPPFAC = 1.0    # potentiation factor
AMPASUPP = 0.4    # fraction of AMPA during storage
STDPTHRESH = -55    # voltage threshold for STDP
STDPSTART = h.STARTDEL+(h.THETA/2)    # STDP starts at same time as EC input
STDPINT = h.THETA/2    # STDP interburst (recall) interval
STDPLEN = h.THETA/2    # STDP burst (storage) length

C_P = 1  # probability of excitatory connections received by each CA1 PC
         # from CA3 inputs (1 gives full connectivity)
def calcSPATT(scaleDown):        
    SPATT = 20*scaleDown	# number of active cells per pattern
    return SPATT
NPATT = 5	# number of stored patterns
NSTORE = 5	# number of new patterns to store

CPATT = 1	# index of cue pattern
CFRAC = 1	# fraction of active cells in cue
iPPC=1		# index of a pattern PC (1st patt in 5 patterns)
iNPPC=0		# index of a non-pattern PC (1st patt in 5 patterns)

# Septal inhibition
SEPNUM = 1000    # number of SEP spikes
SEPSTART = h.STARTDEL+(h.THETA/12)    # time of first SEP spike
SEPINT = 20    # SEP spike ISI (during burst)
SEPNOISE = 0.4    # SEP ISI noise
SEPBINT = 2*h.THETA/3    # SEP interburst interval
SEPBLEN = h.THETA/3    # SEP burst length
SEPWGT = 0.02	# SEP weight to BCs and AACs
SEPWGTL = 0.0002	# SEP weight to BSCs and OLMs
SEPDEL = 1	# SEP delay

# Background excitation (not used)
ENUM = 0    # number of spikes
ESTART = 0    # time of first spike
EINT = 200    # spike ISI
ENOISE = 1    # ISI noise
EWGT = 0.001    # excitatory weights (AMPA)
ENWGT = 0.002    # excitatory weights (NMDA)
EDEL = 1    # delay (msecs)

# EC excitation
ECPATT = 1    # index of output pattern
ECNUM = 1000    # number of EC spikes
ECSTART = h.STARTDEL    # time of first EC spike
ECINT = h.GAMMA    # EC spike ISI
ECNOISE = 0.2    # EC ISI noise
ECWGT = 0.0    # EC weight to PCs
#ECWGT = 0.001    # EC weight to PCs
ECDEL = 1    # EC delay
EIWGT = 0.00015    # excitatory weights to INs
EIDEL = 1    # delay (msecs)

# Cue (CA3) excitation
CNUM = 1000    # number of cue spikes
CSTART = h.STARTDEL+h.ECCA3DEL    # time of first cue spike
CINT = h.GAMMA    # cue spike ISI
CNOISE = 0.2    # cue ISI noise
CHWGT = 0.0015    # cue weight
CLWGT = 0.0005    # unlearnt weight (usually 0)
CNWGT = 0.0005    # excitatory weights (NMDA)
CDEL = 1    # cue delay

# STDP configuration
STDPDFAC = 0    # depression factor
STDPPFAC = 0    # potentiation factor
#STDPDFAC = 0.2    # depression factor
#STDPPFAC = 1.0    # potentiation factor
AMPASUPP = 0.4    # fraction of AMPA during storage
STDPTHRESH = -55    # voltage threshold for STDP
STDPSTART = h.STARTDEL+(h.THETA/2)    # STDP starts at same time as EC input
STDPINT = h.THETA/2    # STDP interburst (recall) interval
STDPLEN = h.THETA/2    # STDP burst (storage) length
