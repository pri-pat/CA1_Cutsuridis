# CA1 heteroassociative memory network: Storage and recall
# CA1 PCs, BCs, AACs, BSCs and OLMs (using moderate cell models)
# EC, CA3 (excitatory) and Septal (inhibitory) inputs
# Cycle is: Recall-Storage-Recall etc
# Serial code adapted from Hines' ran3ser.hoc
# VCU & BPG 8-1-09

# Results reported in V. Cutsuridis, S. Cobb and B.P. Graham,
# "Encoding and retrieval in a model of the hippocampal CA1 microcircuit",
# Hippocampus, in press, DOI 10.1002/hipo.20661, 2009.

from neuron import h, gui
import numpy as np
import random
import math
import sys

h("strdef simname")
h("batchflag = 0")
h("plotflag = 0")
h("scaleDown = 1")
h("scaleEScon = 1")

h.simname="test"
if len(sys.argv)>1:
    h.simname = sys.argv[1]
    if len(sys.argv)>2:
        h.batchflag = int(sys.argv[2])
        if len(sys.argv)>3:
            h.plotflag = int(sys.argv[3])
            if len(sys.argv)>4:
                h.scaleDown = float(sys.argv[4])
                if len(sys.argv)>5:
                    h.scaleEScon = float(sys.argv[5])

print("h.simname = ", h.simname)

h('load_file("ranstream.hoc")')  # to give each cell its own sequence generator

h('strdef fstem')
h.fstem = "Results/" + h.simname
print("simname = ", h.simname, ", fstem = ", h.fstem)


# Set Timing Parameters

h('STARTDEL = 50')	# msecs
h('THETA = 250')	# msecs (4 Hz)
h('GAMMA = 25')	# msecs (40 Hz)
h('ECCA3DEL = 9')	# msecs
h('SIMDUR = 0')	# msecs

# h.SIMDUR = h.STARTDEL + (h.THETA*2) # Scaled down for testing code for technical bugs only
h.SIMDUR = h.STARTDEL + (h.THETA*8)	# simulation duration (msecs)

# Set GID ranges of cells and Load Cell Class definitions
class CellPop:
	def __init__(self, num=0, gidst=0, gidend=0, order=0, filename="", popname="", classtype="", isart=0):
		self.num=num
		self.gidst=gidst
		self.gidend=gidend
		self.order=order
		self.filename=filename
		self.popname=popname
		self.classtype=classtype
		self.isart=isart

poplist=[]
poplist.append(CellPop(num=100*scaleDown, order=0, filename="pyramidal_cell_14Vb.hoc", classtype="PyramidalCell", popname="PyramidalCell", isart=0))
poplist.append(CellPop(num=2, order=1, filename="basket_cell17S.hoc", popname="BasketCell", classtype="BasketCell", isart=0))
poplist.append(CellPop(num=1, order=2, filename="axoaxonic_cell17S.hoc", popname="AACell", classtype="AACell", isart=0))
poplist.append(CellPop(num=1, order=3, filename="bistratified_cell13S.hoc", popname="BistratifiedCell", classtype="BistratifiedCell", isart=0))
poplist.append(CellPop(num=1, order=4, filename="olm_cell2.hoc", popname="OLMCell", classtype="OLMCell", isart=0))
poplist.append(CellPop(num=100*scaleDown, order=5, filename="stim_cell.hoc", popname="CA3Cell", classtype="StimCell", isart=1))
poplist.append(CellPop(num=20*scaleDown, order=6, filename="stim_cell.hoc", popname="ECCell", classtype="StimCell", isart=1))
poplist.append(CellPop(num=10*scaleDown, order=7, filename="burst_cell.hoc", popname="SEPCell", classtype="BurstCell", isart=1))

st=0
for pop in poplist:
	h('load_file("' + pop.filename + '")')
	pop.gidst=st
	pop.gidend=pop.gidst + pop.num - 1
	st = pop.gidend + 1

# Calculate totals of cells
ncell = sum([o.num for o in poplist if o.isart==0])	# total number of real cells
nstim = sum([o.num for o in poplist if o.isart==1])		# total number of inputs
ntot = ncell+nstim # total number of cells

iPC = dictpop['PyramidalCell'].gidst
iBC = dictpop['BasketCell'].gidst
iAAC = dictpop['AACell'].gidst
iBSC = dictpop['BistratifiedCell'].gidst
iOLM = dictpop['OLMCell'].gidst
iCA3 = dictpop['CA3Cell'].gidst
iEC = dictpop['ECCell'].gidst
iSEP = dictpop['SEPCell'].gidst

# Create a dictionary of the cell populations
dictpop = dict(zip([o.popname for o in poplist], poplist))




# creates the cells and appends them to a List called cells
# argument is the number of cells to be created
cells = []
ranlist = []
gidvec = []

i=0
for pop in poplist:
	for j in range(pop.num):	
		exec("cells.append(h." + pop.classtype + ")")
		ranlist.append(h.RandomStream(i))  # ranlist.o(i) corresponds to
		gidvec.append(i)
		i+=1

# Configures the stimulation:


# Septal inhibition
SEPNUM = 1000	# number of SEP spikes
SEPSTART = STARTDEL+(THETA/12)	# time of first SEP spike
SEPINT = 20	# SEP spike ISI (during burst)
SEPNOISE = 0.4	# SEP ISI noise
SEPBINT = 2*THETA/3	# SEP interburst interval
SEPBLEN = THETA/3	# SEP burst length


# Background excitation (not used)
ENUM = 0	# number of spikes
ESTART = 0	# time of first spike
EINT = 200	# spike ISI
ENOISE = 1	# ISI noise
EWGT = 0.001	# excitatory weights (AMPA)
ENWGT = 0.002	# excitatory weights (NMDA)
EDEL = 1	# delay (msecs)

# EC excitation
ECPATT = 1	# index of output pattern
ECNUM = 1000	# number of EC spikes
ECSTART = STARTDEL	# time of first EC spike
ECINT = GAMMA	# EC spike ISI
ECNOISE = 0.2	# EC ISI noise
ECWGT = 0.0	# EC weight to PCs
#ECWGT = 0.001	# EC weight to PCs
ECDEL = 1	# EC delay
EIWGT = 0.00015	# excitatory weights to INs
EIDEL = 1	# delay (msecs)

# Cue (CA3) excitation
CNUM = 1000	# number of cue spikes
CSTART = STARTDEL+ECCA3DEL	# time of first cue spike
CINT = GAMMA	# cue spike ISI
CNOISE = 0.2	# cue ISI noise
CHWGT = 0.0015	# cue weight
CLWGT = 0.0005	# unlearnt weight (usually 0)
CNWGT = 0.0005	# excitatory weights (NMDA)
CDEL = 1	# cue delay

# STDP configuration
STDPDFAC = 0	# depression factor
STDPPFAC = 0	# potentiation factor
#STDPDFAC = 0.2	# depression factor
#STDPPFAC = 1.0	# potentiation factor
AMPASUPP = 0.4	# fraction of AMPA during storage
STDPTHRESH = -55	# voltage threshold for STDP
STDPSTART = STARTDEL+(THETA/2)	# STDP starts at same time as EC input
STDPINT = THETA/2	# STDP interburst (recall) interval
STDPLEN = THETA/2	# STDP burst (storage) length


i=0
for pop in poplist:
	for j in range(pop.num):	
		exec("cells.append(h." + pop.classtype + ")")
		ranlist.append(h.RandomStream(i))  # ranlist.o(i) corresponds to
		gidvec.append(i)
		i+=1

for i in range(dictpop["CA3Cell"].gidst,dictpop["CA3Cell"].gidend+1):
	cells[i].stim.number = ENUM
	cells[i].stim.start = ESTART
	cells[i].stim.interval = EINT
	cells[i].stim.noise = ENOISE

for i in range(dictpop["ECCell"].gidst,dictpop["ECCell"].gidend+1):
	cells[i].stim.number = ENUM
	cells[i].stim.start = ESTART
	cells[i].stim.interval = EINT
	cells[i].stim.noise = ENOISE

for i in range(dictpop["SEPCell"].gidst,dictpop["SEPCell"].gidend+1):
	cells[i].stim.number = SEPNUM
	cells[i].stim.start = SEPSTART
	cells[i].stim.interval = SEPINT
	cells[i].stim.noise = SEPNOISE
	cells[i].stim.burstint = SEPBINT
	cells[i].stim.burstlen = SEPBLEN

    # Use the gid-specific random generator so random streams are
    # independent of where and how many stims there are.
    h('rs = ranlist[i]')
    h('stim.noiseFromRandom(rs.r)')
    h('rs.r.negexp(1)')
    h('rs.start()')


# Pattern storage and recall parameters

C_P = 1  # probability of excitatory connections received by each CA1 PC
         # from CA3 inputs (1 gives full connectivity)
         
SPATT = 20*scaleDown	# number of active cells per pattern
NPATT = 5	# number of stored patterns
NSTORE = 5	# number of new patterns to store

CPATT = 1	# index of cue pattern
CFRAC = 1	# fraction of active cells in cue
iPPC=1		# index of a pattern PC (1st patt in 5 patterns)
iNPPC=0		# index of a non-pattern PC (1st patt in 5 patterns)

# file name of connection weights and patterns
# (cue and EC patterns taken from FSTORE file to implement storage)
# (use same file for FPATT and FSTORE to test recall only)
if scaleDown==1:
	FCONN = "weights/wgtsN100S20P5.dat" #"Weights/wgtsN100S20P5.dat"
	FPATT = "Weights/pattsN100S20P5.dat"	# "Weights/pattsN100S20P5.dat"	# already stored patterns
	FSTORE = "Weights/pattsN100S20P5.dat"	# "Weights/pattsN100S20P5.dat"	# new patterns to store
else:
	FCONN = "weights/wgtsN100S20P5scaled.dat" #"Weights/wgtsN100S20P5.dat"
	FPATT = "Weights/pattsN100S20P5scaled.dat"	# "Weights/pattsN100S20P5.dat"	# already stored patterns
	FSTORE = "Weights/pattsN100S20P5scaled.dat"	# "Weights/pattsN100S20P5.dat"	# new patterns to store

# Connectivity
# of connections received by each POSTCELLTYPE from a given PRECELLTYPE
class popConn:
	def __init__(self, popname="", prepop="", prenum=1, type="", weight=0, delay=3, synst=0, synend=0):
		self.popname=popname
		self.prepop=prepop
		self.prenum=prenum
		self.type=type
		self.weight=weight
		self.delay=delay
		self.synst=synst
		self.synend=synend

# Synapse indices
# onto CA1 PCs
E_EC = 0	# EC AMPA excit to medium SLM (2 of)
E_CA3 = 2	# CA3 AMPA excit to medium SR
EN_CA3 = 3	# CA3 NMDA excit to medium SR
EM_CA3 = 23	# CA3 modifiable (STDP) AMPA excit to medium SR
E_PC = 4	# CA1 recurrent AMPA excit to proximal SR
I_BC = 5	# ff&fb inhib via BCs to soma
I_AAC = 6	# ff&fb inhib via AACs to axon initial segment
I_BSC = 11	# ff&fb inhib via BSCs to SR med (12 of: 6 GABAA, 6 GABAB)
I_OLM = 7	# fb inhib via OLMs to SLM (4 of: 2 GABAA, 2 GABAB)

# onto INs (BC, AAC, BSC)
EI_EC = 0	# EC AMPA excit (2 of; not onto BSC)
EI_CA3 = 2	# CA3 AMPA excit (4 of)
EI_PC = 6	# CA1 PC AMPA excit (2 of)
II_SAME = 8	# inhib from neighbouring INs (BC->BC; BSC->BSC)
II_OPP = 9	# inhib from other INs (BSC->BC; BC->BSC)
II_SEP = 10	# inhib from septum (4 of: 2 GABAA, 2 GABAB)

# onto INs (OLM)
EO_PC = 0	# CA1 PC AMPA excit (2 of)
IO_IN = 2	# inhib from INs and septum (2 of: 1 GABAA, 1 GABAB)

connlist=[]                    # Postsynaptic type      Presynaptic type  number of connections (at least 1)
connlist.append(popConn(popname="PyramidalCell", prepop="CA3Cell", prenum=max([nCA3 * scaleEScon, 1]), type="AMPA", weight=0, delay=1, synst=, synend=)) # CA3_PC
connlist.append(popConn(popname="BasketCell",    prepop="CA3Cell", prenum=max([nCA3 * scaleEScon, 1]), type="AMPA")) # CA3_BC
connlist.append(popConn(popname="AACell",        prepop="CA3Cell", prenum=max([nCA3 * scaleEScon, 1]), type="AMPA")) # CA3_AAC
connlist.append(popConn(popname="BistratifiedCell", prepop="CA3Cell", prenum=max([nCA3 * scaleEScon, 1]), type="AMPA")) # CA3_BSC

connlist.append(popConn(popname="PyramidalCell", prepop="ECCell", prenum=max([nEC * scaleEScon, 1]), type="AMPA")) # EC_PC
connlist.append(popConn(popname="BasketCell",    prepop="ECCell", prenum=max([nEC * scaleEScon, 1]), type="AMPA")) # EC_BC
connlist.append(popConn(popname="AACell",        prepop="ECCell", prenum=max([nEC * scaleEScon, 1]), type="AMPA")) # EC_AAC

connlist.append(popConn(popname="BasketCell",    prepop="SEPCell", prenum=nSEP, type="GABAA", weight=0.02, delay=1)) # SEP_BC
connlist.append(popConn(popname="AACell",        prepop="SEPCell", prenum=nSEP, type="GABAA", weight=0.02, delay=1)) # SEP_AAC
connlist.append(popConn(popname="BistratifiedCell", prepop="SEPCell", prenum=nSEP, type="GABAA", weight=0.0002, delay=1)) # SEP_BSC
connlist.append(popConn(popname="OLMCell",       prepop="SEPCell", prenum=nSEP, type="GABAA", weight=0.0002, delay=1)) # SEP_OLM

connlist.append(popConn(popname="PyramidalCell", prepop="PyramidalCell", prenum=max([1 * scaleEScon, 1]), type="AMPA", weight=0.001, delay=1)) # PC_PC
connlist.append(popConn(popname="BasketCell",    prepop="PyramidalCell", prenum=max([dictpop['PyramidalCell'].num * scaleEScon, 1]), type="AMPA", weight = 0.0005, delay = 1)) # PC_BC
connlist.append(popConn(popname="BistratifiedCell", prepop="PyramidalCell", prenum=max([dictpop['PyramidalCell'].num * scaleEScon, 1]), type="AMPA", weight = 0.0005, delay = 1)) # PC_BSC
connlist.append(popConn(popname="AACell",        prepop="PyramidalCell", prenum=max([dictpop['PyramidalCell'].num * scaleEScon, 1]), type="AMPA", weight = 0.0005, delay = 1)) # PC_AAC
connlist.append(popConn(popname="OLMCell",       prepop="PyramidalCell", prenum=max([dictpop['PyramidalCell'].num * scaleEScon, 1]), type="AMPA", weight = 0.00005, delay = 1)) # PC_OLM

connlist.append(popConn(popname="PyramidalCell", prepop="BasketCell", prenum=2, type="GABAA", weight = 0.02, delay = 1)) # BC_PC
connlist.append(popConn(popname="BistratifiedCell", prepop="BasketCell", prenum=2, type="GABAA", weight = 0.02, delay = 1)) # BC_BSC

# The call to make this connection is commented out in hoc code
# connlist.append(popConn(popname="OLMCell",       prepop="BasketCell", prenum=2, type="GABAA", weight = 0.0, delay = 1)) # BC_OLM
connlist.append(popConn(popname="BasketCell",    prepop="BasketCell", prenum=1, type="GABAA", weight = 0.001, delay = 1)) # BC_BC

connlist.append(popConn(popname="PyramidalCell", prepop="AACell", prenum=1, type="GABAA", weight = 0.04, delay = 1)) # AAC_PC
connlist.append(popConn(popname="PyramidalCell", prepop="BistratifiedCell", prenum=1, type="GABAA", weight = 0.002, delay = 1)) # BSC_PC
connlist.append(popConn(popname="PyramidalCell", prepop="BistratifiedCell", prenum=1, type="GABAB", weight = 0.0004, delay = 1)) # BSC_PC

connlist.append(popConn(popname="BasketCell",    prepop="BistratifiedCell", prenum=1, type="GABAA", weight = 0.01, delay = 1)) # BSC_BC

connlist.append(popConn(popname="PyramidalCell", prepop="OLMCell", prenum=1, type="GABAA", weight = 0.04, delay = 1)) # OLM_PC
connlist.append(popConn(popname="PyramidalCell", prepop="OLMCell", prenum=1, type="GABAB", weight = 0.0004, delay = 1)) # OLM_PC

# In the original Cutsuridis code, this weight is set first to 0.01 and then overwritten a few lines later to 0.0.
# The call to make this conection is commented out in hoc
# connlist.append(popConn(popname="BasketCell",    prepop="OLMCell", prenum=1, type="GABAA", weight = 0.0, delay = 1)) # OLM_BC
  
connect_random_low_start_ = 1  # low seed for mcell_ran4_init()
h.mcell_ran4_init(connect_random_low_start_)
nclist = []

# Make connections with data from above
for conn in connlist: 
	connectcells(dictpop[conn.popname].num, dictpop[conn.popname].gidst, dictpop[conn.prepop].num, dictpop[conn.prepop].gidst, conn.prenum, conn.synst, conn.synend, conn.delay, conn.weight)


# print "   EC..."
# EC to BC
connectcells(nbcell, iBC, nEC, iEC, EC_BC, EI_EC, EI_EC+1, EIDEL, EIWGT)
# EC to AAC
connectcells(naacell, iAAC, nEC, iEC, EC_AAC, EI_EC, EI_EC+1, EIDEL, EIWGT)

#print "   CA3..."
# CA3 to BC
connectcells(nbcell, iBC, nCA3, iCA3, CA3_BC, EI_CA3, EI_CA3+3, EIDEL, EIWGT)
# CA3 to AAC
connectcells(naacell, iAAC, nCA3, iCA3, CA3_AAC, EI_CA3, EI_CA3+3, EIDEL, EIWGT)
# CA3 to BSC
connectcells(nbscell, iBSC, nCA3, iCA3, CA3_BSC, EI_CA3, EI_CA3+3, EIDEL, EIWGT)

#print "   septum..."
# SEP to BC
connectcells(nbcell, iBC, nSEP, iSEP, SEP_BC, II_SEP, II_SEP+1, SEPDEL, SEPWGT)
# SEP to AAC
connectcells(naacell, iAAC, nSEP, iSEP, SEP_AAC, II_SEP, II_SEP+1, SEPDEL, SEPWGT)
# SEP to BSC
connectcells(nbcell, iBC, nSEP, iSEP, SEP_BSC, II_SEP, II_SEP+1, SEPDEL, SEPWGTL)
# SEP to OLM
connectcells(naacell, iAAC, nSEP, iSEP, SEP_OLM, IO_IN, IO_IN, SEPDEL, SEPWGTL)

#print "   PC..."
# PC to PC
connectcells(npcell, iPC, npcell, iPC, PC_PC, E_PC, E_PC, Pcell2Pcell_delay, Pcell2Pcell_weight)
# PC to BC
connectcells(nbcell, iBC, npcell, iPC, PC_BC, EI_PC, EI_PC+1, Pcell2Bcell_delay, Pcell2Bcell_weight)
# PC to AAC
connectcells(naacell, iAAC, npcell, iPC, PC_AAC, EI_PC, EI_PC+1, Pcell2AAcell_delay, Pcell2AAcell_weight)
# PC to BSC
connectcells(nbscell, iBSC, npcell, iPC, PC_BSC, EI_PC, EI_PC+1, Pcell2BScell_delay, Pcell2BScell_weight)
# PC to OLM
connectcells(nolm, iOLM, npcell, iPC, PC_OLM, EO_PC, EO_PC+1, Pcell2OLMcell_delay, Pcell2OLMcell_weight)

#print "   INs..."
# BC to PC
connectcells(npcell, iPC, nbcell, iBC, BC_PC, I_BC, I_BC, Bcell2Pcell_delay, Bcell2Pcell_weight)
# BC to BC
connectcells(nbcell, iBC, nbcell, iBC, BC_BC, II_SAME, II_SAME, Bcell2Bcell_delay, Bcell2Bcell_weight)
# BC to BSC
connectcells(nbscell, iBSC, nbcell, iBC, BC_BSC, II_OPP, II_OPP, Bcell2BScell_delay, Bcell2BScell_weight)
# BC to OLM - unsure why commented out when defined above
#connectcells(nolm, iOLM, nbcell, iBC, BC_OLM, IO_IN, IO_IN, Bcell2OLMcell_delay, Bcell2OLMcell_weight)
# AAC to PC
connectcells(npcell, iPC, naacell, iAAC, AAC_PC, I_AAC, I_AAC, AAcell2Pcell_delay, AAcell2Pcell_weight)
# BSC to PC
connectcells(npcell, iPC, nbscell, iBSC, BSC_PC, I_BSC, I_BSC+5, BScell2Pcell_delay, BScell2Pcell_weight)
connectcells(npcell, iPC, nbscell, iBSC, BSC_PC, I_BSC+6, I_BSC+11, BScell2Pcell_delay, BScell2Pcell_GABAB_weight)
# BSC to BSC - not defined above and commented out in hoc code
#connectcells(nbscell, iBSC, nbscell, iBSC, BSC_BSC, II_SAME, II_SAME, BScell2BScell_delay, BScell2BScell_weight)
# BSC to BC
connectcells(nbcell, iBC, nbscell, iBSC, BSC_PC, II_OPP, II_OPP, BScell2Bcell_delay, BScell2Bcell_weight)
# OLM to PC
connectcells(npcell, iPC, nolm, iOLM, OLM_PC, I_OLM, I_OLM+1, OLMcell2Pcell_delay, OLMcell2Pcell_weight)
connectcells(npcell, iPC, nolm, iOLM, OLM_PC, I_OLM+2, I_OLM+3, OLMcell2Pcell_delay, OLMcell2Pcell_GABAB_weight)
# OLM to BC - commented out, defined above and then weight set to 0 later
#connectcells(nbcell, iBC, nolm, iOLM, OLM_BC, II_OPP, II_OPP, OLMcell2Bcell_delay, OLMcell2Bcell_weight)




h.xopen("HAM_StoRec_ser_diet.hoc")