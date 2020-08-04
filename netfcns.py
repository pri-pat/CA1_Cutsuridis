# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 00:37:03 2020

@author: mbezaire
"""
import numpy as np
import random
from neuron import h
from model_const import *
import matplotlib.pyplot as plt

def connectcells(cells, ranlist, nclist, npost, postgidstart,npre,pregidstart,synstart,synend,npresyn,weight,delay): # {local i, j, gid, nsyn, r  localobj syn, nc, rs, u
    # initialize the pseudorandom number generator
    
    gids = [x.gid for x in cells]
    
    for gid in range(int(postgidstart),int(postgidstart)+int(npost)): # appropriate target cell
        i=gids.index(gid)
        rs = ranlist[i]  # RandomStream for cells.object(i)
        rs.start()
        rs.r.discunif(int(pregidstart),int(pregidstart)+npre-1)  # pick a random presynaptic cell by gid of presynaptic cell type return source cell index
        u = np.zeros(int(npre))  # for sampling without replacement, u[i]==1 means spike source i has already been chosen
        nsyn = 0
        
        cell = cells[i]
        while (nsyn < npresyn and nsyn < npre):
            r = int(rs.repick())
            # no self-connection and only one connection from any source
            if (r != cell.gid and u[r-int(pregidstart)] == 0):
                # target synapses
                for j in range(synstart,synend+1) :
                    # set up connection from source to target
                    syn = cells[i].pre_list[j]
                    #nc = pc.gid_connect(r, syn)
                    nc = cells[r].connect2target(syn)
                    nclist.append(nc)
                    nc.delay = delay
                    nc.weight[0] = weight
                
                u[r-int(pregidstart)] = 1
                nsyn += 1


# connects the EC input layer to PC cells
# read active PCs from pattern file
# all-to-all connectivity between EC and PC pattern cells
# appends the PC NetCons to a List called ncslist
def connectEC(FPATT, ECPATT, NPATT, synstart, numsyn, cells, iPC, iEC, dictpop):# {local i, gid, ncue  localobj cue, cstim, syn, src, nc, fp, target
    ncelist = []
    # read pattern file (ECPATT=num rows, NPATT = num columns)
    cue = np.loadtxt(fname = FPATT)
    if (cue.shape != (ECPATT, NPATT)):
        print("The cue data is a different shape than expected:", cue.shape)

    # find active cells in pattern
    for i in range(len(cue)):
        #if (!pc.gid_exists(i+iPC)) { continue }
        if (cue[i,0] == 1): # TODO added a column index that wouldn't be there usually?
            # print "Pattern cell ", i
            # target = pc.gid2cell(i+iPC)
            target = cells[i+iPC]
            #target synapses
            for k in range(synstart, synstart+numsyn):
                syn = target.pre_list[k]    # excitatory synapse
                # create pattern stimulus
                for j in range(int(dictpop["ECCell"].num)):
                    src = cells[int(j+iEC)].stim
                    # set up connection from source to target
                    #nc = pc.gid_connect(j+iEC, syn)
                    nc = h.NetCon(src, syn)
                    ncelist.append(nc)
                    nc.delay = ECDEL
                    nc.weight[0] = ECWGT
                    

# connects the CA3 input layer to output cells (PCs and INs)
# read PC connections from a file, with connections to
# a target being a column with index i for target cell i
# appends the PC NetCons to a List called ncslist
                    
def connectCA3(FCONN, C_P, EM_CA3, EN_CA3, cells, dictpop): # {local i, j, cp, gid  localobj src, syn, synN, nc, fc, rs, conns, rc
    cp = C_P    # connection probability
    #mcell_ran4_init(connect_random_low_start_)
    #rc = new Vector(dictpop["CA3Cell"].num)  # random physical connectivity
    ncslist = []
    random.seed(connect_random_low_start_)
    #ncilist = new List()
    # inputs to PCs determined by weight matrix
    for cell in cells:    # loop over possible target cells
        gid = cell.gid    # id of cell
        if (gid >= dictpop['PyramidalCell'].gidst and gid < int(dictpop["PyramidalCell"].num)+dictpop['PyramidalCell'].gidst):    # appropriate target cell
            syn = cell.pre_list[EM_CA3]    # AMPA synapse with STDP
            syn.wmax = CHWGT
            syn.wmin = CLWGT
            syn.d = STDPDFAC    # depression
            syn.p = STDPPFAC    # potentiation
            syn.gscale = AMPASUPP    # fraction of AMPA during storage
            syn.thresh = STDPTHRESH    # threshold for postsynaptic voltage detection
            syn.gbdel = STDPSTART
            syn.gbint = STDPINT
            syn.gblen = STDPLEN
            synN = cell.pre_list[EN_CA3]    # NMDA synapse
            #rs = ranlist[i]  # the corresponding RandomStream
            #rs.start()
            #rs.r.uniform(0, 1)  # return integer in range 0..1
            #rc.setrand(rs.r)    # generate random connectivity
            
            # open connections file
            # read incoming weights for cell gid
            conns = np.loadtxt(fname = FCONN)
    
        for j in range(int(dictpop["CA3Cell"].num)):
            #for j=0, CA3_PC-1 { # You might need something like this line instead so that it doesn't error out if you are only planning to make a scaled down number of synapses
            # only connection if physical connection exists
            #if (rc[j] <= cp):
            if (random.uniform(0,1) <= cp):
                #print "   src ", j
                src = cells[int(j+dictpop['CA3Cell'].gidst)].stim
                # set up connection from source to target NMDA synapse
                #        nc = pc.gid_connect(j+iCA3, synN)
                nc = h.NetCon(src, synN)
                ncslist.append(nc)
                nc.delay = CDEL
                nc.weight[0] = CNWGT    # NMDA weight same for all connections
                # high AMPA if weight is 1
                if (conns[j,0] == 1): # TODO access the correct column
                    # set up connection from source to target
                    #nc = pc.gid_connect(j+iCA3, syn)
                    nc = h.NetCon(src, syn)
                    ncslist.append(nc)
                    nc.delay = CDEL
                    nc.weight[0] = CHWGT
                else:
                    # set up connection from source to target
                    #nc = pc.gid_connect(j+iCA3, syn)
                    nc = h.NetCon(src, syn)
                    ncslist.append(nc)
                    nc.delay = CDEL
                    nc.weight[0] = CLWGT    # unlearned weight


# sets the CA3, EC and Septal background inputs
def mkinputs(cells, iCA3, iEC, iSEP, ntot, dictpop): #{local i localobj stim, rs
    for cell in cells:
        gid = cell.gid    # id of cell
        if (gid >= iCA3 and gid < ntot-dictpop["SEPCell"].num-dictpop["ECCell"].num):    # appropriate target cell
            # set background activity for excitatory inputs
            stim = cell.stim
            stim.number = ENUM
            stim.start = ESTART
            stim.interval = EINT
            stim.noise = ENOISE
        
        if (gid >= iEC and gid < ntot-dictpop["SEPCell"].num):    # appropriate target cell
            # set background activity for excitatory inputs
            stim = cell.stim
            stim.number = ENUM
            stim.start = ESTART
            stim.interval = EINT
            stim.noise = ENOISE
    
        if (gid >= iSEP and gid < ntot):    # appropriate target cell
            # set background activity for septum
            stim = cell.stim
            rs = ranlist[i]
            stim.number = SEPNUM
            stim.start = SEPSTART
            stim.interval = SEPINT
            stim.noise = SEPNOISE
            stim.burstint = SEPBINT
            stim.burstlen = SEPBLEN
            # Use the gid-specific random generator so random streams are
            # independent of where and how many stims there are.
            stim.noiseFromRandom(rs.r)
            rs.r.negexp(1)
            rs.start()
            
#########################
# Instrumentation, i.e. stimulation and recording
#########################

# setup activity in EC stims
def mkEC(cells, ranlist, iEC, nEC): # {local i, necs localobj cstim, rs
    EClist = []
    necs = 0
    print("Make EC input...")
    for i, cell in enumerate(cells):
        gid = cell.gid    # id of cell
        if (gid >= iEC and gid < iEC+nEC):     # appropriate target cell
            # create cue stimulus
            cstim = cell.stim
            rs = ranlist[i]
            cstim.number = ECNUM
            cstim.start = ECSTART
            cstim.interval = ECINT
            cstim.noise = ECNOISE
            # Use the gid-specific random generator so random streams are
            # independent of where and how many stims there are.
            cstim.noiseFromRandom(rs.r)
            rs.r.normal(0, 1)
            rs.start()
            EClist.append(i)
            necs += 1




# setup activity pattern in input cue stims
def mkcue(FPATT, CPATT, CFRAC, NPATT, SPATT, cells, iCA3, ranlist):
    print( "Make cue (CA3) input...")
    cuelist = []
    # open patterns file
    cue = np.loadtxt(fname = FPATT) # read pattern


    ncue = 0
    # find active cells in pattern
    for i in range(len(cue)):
        #if (!pc.gid_exists(i+iCA3)) { continue }
        if (ncue <= SPATT*CFRAC):     # fraction of active cells in cue
            if (cue[i,0] == 1): #TODO find the correct column
                print("Cue cell ", i)
                #cstim = pc.gid2cell(i+iCA3)
                cstim = cells[int(i+iCA3)].stim
                for j, cell in enumerate(cells):
                    if (cell.gid == i+iCA3):
                        break    # find cell index
                    
                    rs = ranlist[j]
                    # create cue stimulus
                    cstim.number = CNUM
                    cstim.start = CSTART
                    cstim.interval = CINT
                    cstim.noise = CNOISE
                    # Use the gid-specific random generator so random streams are
                    # independent of where and how many stims there are.
                    cstim.noiseFromRandom(rs.r)
                    rs.r.normal(0, 1)
                    rs.start()
                    cuelist.append(i)
                    ncue += 1
                    
                    print("  cue size ", ncue)

# remove activity pattern in input cue stims
def erasecue(): # {local i, j localobj cstim
  for i in range(len(cuelist)):
    #if (!pc.gid_exists(i+iCA3)) { continue }
    #cstim = pc.gid2cell(i+iCA3)
    cstim = cells[cuelist[i]+iCA3].stim
    cstim.number = 0

# Spike recording
# tvec, idvec will be Vectors that record all spike times (tvec)
# and the corresponding id numbers of the cells that spiked (idvec)
def spikerecord(cells):
    print( "Record spikes...")
    for cell in cells:
        if (cell.is_art==0):
            cell._spike_detector = h.NetCon(cell.soma(0.5)._ref_v, None, sec=cell.soma)
            cell.spike_times = h.Vector()
            cell._spike_detector.record(cell.spike_times)
            cell.soma_v = h.Vector().record(cell.soma(0.5)._ref_v)
        #else: # to print stim cell spikes...


# Record cell voltage traces
# Vectors that record voltages from pattern PC
# Vectors that record voltages from non-pattern PC
# Vectors that record voltages from INs

def vrecord(cells,dictpop):
    iPC = dictpop['PyramidalCell'].gidst
    iBC = dictpop['BasketCell'].gidst
    iAAC = dictpop['AACell'].gidst
    iBSC = dictpop['BistratifiedCell'].gidst
    iOLM = dictpop['OLMCell'].gidst
    iCA3 = dictpop['CA3Cell'].gidst
    iEC = dictpop['ECCell'].gidst
    iSEP = dictpop['SEPCell'].gidst    
    
    print( "Record example voltage traces...")
    results = {}
    for cell in cells:	# loop over possible target cells
        gid = cell.gid	# id of cell
        if (gid==iPPC):
            results["pvsoma"] = h.Vector().record(cell.soma(0.5)._ref_v)
            results["pvsr"] = h.Vector().record(cell.radTmed(0.5)._ref_v)
            results["pvslm"] = h.Vector().record(cell.lm_thick1(0.5)._ref_v)

        if (gid==iNPPC):
            results["npvsoma"] = h.Vector().record(cell.soma(0.5)._ref_v)
            results["npvsr"] = h.Vector().record(cell.radTmed(0.5)._ref_v)
            results["npvslm"] = h.Vector().record(cell.lm_thick1(0.5)._ref_v)

        if (gid==iBC):
            results["vBC"] = h.Vector().record(cell.soma(0.5)._ref_v)
            print("Recording results into vBC from ", cell)

        if (gid==iAAC):
            results["vAAC"] = h.Vector().record(cell.soma(0.5)._ref_v)
            print("Recording results into vAAC from ", cell)

        if (gid==iBSC):
            results["vBSC"] = h.Vector().record(cell.soma(0.5)._ref_v)
            print("Recording results into vBSC from ", cell)

        if (gid==iOLM):
            results["vOLM"] = h.Vector().record(cell.soma(0.5)._ref_v)
            print("Recording results into vOLM from ", cell)

    return results

def spikeout(cells,fstem):
    #printf("\ntime\t cell\n")  # print header once    
    with open("{}_spt.dat".format(fstem), 'w') as f:
        for i, cell in enumerate(cells):
            if (cell.is_art==0):
                for spk in cell.spike_times:
                    f.write("{}\t{}\n".format(spk, cell.gid))

def vout(cells,results,fstem,dictpop):
    

    with open("{}_pvsoma.dat".format(fstem), 'w') as f:
        for v in results["pvsoma"]:
            f.write("{}\n".format(v))

    with open("{}_pvsr.dat".format(fstem), 'w') as f:
        for v in results["pvsr"]:
            f.write("{}\n".format(v))

    with open("{}_pvslm.dat".format(fstem), 'w') as f:
        for v in results["pvslm"]:
            f.write("{}\n".format(v))

    with open("{}_npvsoma.dat".format(fstem), 'w') as f:
        for v in results["npvsoma"]:
            f.write("{}\n".format(v))

    with open("{}_npvsr.dat".format(fstem), 'w') as f:
        for v in results["npvsr"]:
            f.write("{}\n".format(v))

    with open("{}_npvslm.dat".format(fstem), 'w') as f:
        for v in results["npvslm"]:
            f.write("{}\n".format(v))

    with open("{}_BC.dat".format(fstem), 'w') as f:
        for v in results["vBC"]:
            f.write("{}\n".format(v))

    with open("{}_AAC.dat".format(fstem), 'w') as f:
        for v in results["vAAC"]:
            f.write("{}\n".format(v))
        plt.figure()
        plt.plot(results["vAAC"])
        plt.xlabel("AAC")
        plt.show()

    with open("{}_BSC.dat".format(fstem), 'w') as f:
        for v in results["vBSC"]:
            f.write("{}\n".format(v))

    with open("{}_OLM.dat".format(fstem), 'w') as f:
        for v in results["vOLM"]:
            f.write("{}\n".format(v))


# produce raster plot of spiking activity
def spikeplot(cells,tstop,ntot):
    plt.figure()
    for i, cell in enumerate(cells):
        if (cell.is_art==0 and cell.spike_times.size()>0):
            plt.vlines(cell.spike_times, i + 0.5, i + 1.5)
    plt.xlabel('Time (ms)')
    plt.ylabel('Neuron (gid)')
    plt.show()
    
def vplot(cells,results):
    plt.figure()
    t = np.arange(0,h.t,h.dt)
    for key in results:
        plt.plot(t,results[key],label=key)
    plt.xlabel('Time (ms)')
    plt.ylabel('Membrane Potential (mV)')
    plt.legend()
    plt.show()

# panel for simulation results
def xspikeres():
    h.xpanel("Spike results")
    h.xbutton("Write out", "spikeout()")
    h.xbutton("Plot", "spikeplot()")
    h.xpanel()



