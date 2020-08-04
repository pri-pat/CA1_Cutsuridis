# Data from Saraga et al. (2003) paper
# changed morphology and some channel densities (BPG 12-1-09)
#   OLM_Cell
# basic_shape,
# pre_list, connect2target
# soma, dend1, dend2, axon
# all
from neuron import h

class BistratifiedCell():
    """ Bistratified Cell definition """
    def __init__(self, gid = -1):
        self.x = 0; self.y = 0; self.z = 0
        self.gid = gid
        self.create_sections() 
        self.build_topology()
        self.build_subsets() # subsets()
        self.define_geometry() # geom()
        self.define_biophysics() # biophys()
        # pre_list = new List()
        self.addSynapses() # synapses
        self.is_art = 0
        self.nc = []
        
    def __repr__(self):
        return "Bistratified Cell {}".format(self.gid)

    def create_sections(self):
        self.soma = h.Section(name='soma', cell=self)
        self.radT2 = h.Section(name='radT2', cell=self)
        self.radM2 = h.Section(name='radM2', cell=self)
        self.radt2 = h.Section(name='radt2', cell=self)
        self.radT1 = h.Section(name='radT1', cell=self)
        self.radM1 = h.Section(name='radM1', cell=self)
        self.radt1 = h.Section(name='radt1', cell=self)
        self.oriT1 = h.Section(name='oriT1', cell=self)
        self.oriM1 = h.Section(name='oriM1', cell=self)
        self.orit1 = h.Section(name='orit1', cell=self)
        self.oriT2 = h.Section(name='oriT2', cell=self)
        self.oriM2 = h.Section(name='oriM2', cell=self)
        self.orit2 = h.Section(name='orit2', cell=self)

    def build_topology(self):
        self.radT2.connect(self.soma(1))
        self.radM2.connect(self.radT2(1))
        self.radt2.connect(self.radM2(1))
        self.radT1.connect(self.soma(0))
        self.radM1.connect(self.radT1(1))
        self.radt1.connect(self.radM1(1))

        self.oriT1.connect(self.soma(0))
        self.oriM1.connect(self.oriT1(1))
        self.orit1.connect(self.oriM1(1))

        self.oriT2.connect(self.soma(1))
        self.oriM2.connect(self.oriT2(1))
        self.orit2.connect(self.oriM2(1))



    def define_geometry(self):
        for sec in self.all:
            sec.pt3dclear()

        h.pt3dadd(0, 0, 0, 10, sec=self.soma)
        h.pt3dadd(15, 0, 0, 10, sec=self.soma)
        
        h.pt3dadd(15, 0, 0, 1, sec=self.radT2)
        h.pt3dadd(45, 30, 0, 1, sec=self.radT2)
        
        h.pt3dadd(45, 30, 0, 1, sec=self.radM2)
        h.pt3dadd(75, 60, 0, 1, sec=self.radM2)
        
        h.pt3dadd(75, 60, 0, 1, sec=self.radt2)
        h.pt3dadd(90, 75, 0, 1, sec=self.radt2)
        
        h.pt3dadd(0, 0, 0, 1, sec=self.radT1)
        h.pt3dadd(-14, 15, 0, 1, sec=self.radT1)

        h.pt3dadd(-14, 15, 0, 1, sec=self.radM1)
        h.pt3dadd(-29, 30, 0, 1, sec=self.radM1)

        h.pt3dadd(-29, 30, 0, 1, sec=self.radt1)
        h.pt3dadd(-44, 45, 0, 1, sec=self.radt1)

        h.pt3dadd(0, 0, 0, 1, sec=self.oriT1)
        h.pt3dadd(-29, -29, 0, 1, sec=self.oriT1)

        h.pt3dadd(-29, -29, 0, 1, sec=self.oriM1)
        h.pt3dadd(-59, -59, 0, 1, sec=self.oriM1)

        h.pt3dadd(-59, -59, 0, 1, sec=self.orit1)
        h.pt3dadd(-89, -89, 0, 1, sec=self.orit1)

        h.pt3dadd(15, 0, 0, 1, sec=self.oriT2)
        h.pt3dadd(45, -29, 0, 1, sec=self.oriT2)

        h.pt3dadd(45, -29, 0, 1, sec=self.oriM2)
        h.pt3dadd(75, -59, 0, 1, sec=self.oriM2)

        h.pt3dadd(75, -59, 0, 1, sec=self.orit2)
        h.pt3dadd(105, -89, 0, 1, sec=self.orit2)
       
        
        self.soma.L = 20
        self.soma.diam = 10
        self.radT2.L = 100
        self.radT2.diam = 4
        self.radM2.L = 100
        self.radM2.diam = 3
        self.radt2.L = 200
        self.radt2.diam = 2

        self.radT1.L = 100
        self.radT1.diam = 4
        self.radM1.L = 100
        self.radM1.diam = 3
        self.radt1.L = 200
        self.radt1.diam = 2

        self.oriT1.L = 100
        self.oriT1.diam = 2
        self.oriM1.L = 100
        self.oriM1.diam = 1.5
        self.orit1.L = 100
        self.orit1.diam = 1

        self.oriT2.L = 100
        self.oriT2.diam = 2
        self.oriM2.L = 100
        self.oriM2.diam = 1.5
        self.orit2.L = 100
        self.orit2.diam = 1

        
        for sec in self.all:
            h("lf = lambda_f(100)")
            sec.nseg = int((sec.L/(0.1*h.lf)+.9)/2)*2 + 1 

    def build_subsets(self):
        self.all = h.SectionList()
        self.all.wholetree(sec=self.soma)


    def define_biophysics(self):

        gna = 0.3
        
        self.soma.insert("ichan2")
        for seg in self.soma:
            seg.gnatbar_ichan2 = gna  		# 0.12 //original 0.030 to .055 
            seg.gkfbar_ichan2 = 0.013  		# original 0.015
            seg.gl_ichan2 = 0.00018
            cm=1.4
        
        self.radt1.insert("ichan2")
        for seg in self.radt1:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
            cm=1.4

        
        self.radt2.insert("ichan2")
        for seg in self.radt2:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
            cm=1.4

        
        self.radM1.insert("ichan2")
        for seg in self.radM1:
            seg.gnatbar_ichan2 = gna  		# 0.3  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
            cm=1.4
        
        self.radM2.insert("ichan2")
        for seg in self.radM2:
            seg.gnatbar_ichan2 = gna  		# 0.3  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
            cm=1.4

        
        self.radT1.insert("ichan2")
        for seg in self.radT1:
            seg.gnatbar_ichan2 = gna  		# 0.2  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
            cm=1.4
        
        self.radT2.insert("ichan2")
        for seg in self.radT2:
            seg.gnatbar_ichan2 = gna  		# 0.2  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
            cm=1.4
        
        self.oriT1.insert("ichan2")
        for seg in self.oriT1:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
            cm=1.4
        
        self.oriT2.insert("ichan2")
        for seg in self.oriT2:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
            cm=1.4

        
        self.oriM1.insert("ichan2")
        for seg in self.oriM1:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
            cm=1.4
        
        self.oriM2.insert("ichan2")
        for seg in self.oriM2:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
            cm=1.4
        
        self.orit1.insert("ichan2")
        for seg in self.orit1:
            seg.gnatbar_ichan2 = gna  		# Sodium conductance (original 0.015)
            seg.gkfbar_ichan2 = 0.013  		# Delayed K+ rectifier (fast)
            seg.gl_ichan2 = 0.00018         # Leak conductance
            cm=1.4
        
        self.orit2.insert("ichan2")
        for seg in self.orit2:
            seg.gnatbar_ichan2 = gna  		# Sodium conductance (original 0.015)
            seg.gkfbar_ichan2 = 0.013  		# Delayed K+ rectifier (fast)
            seg.gl_ichan2 = 0.00018         # Leak conductance
            cm=1.4
        
        for sec in self.all:
            # self.cm = Not setting cm
            self.Ra = 100			# 31.3 +/- 10.9
            self.enat = 55
            self.ekf = -90
            self.ek = -90
            self.elca = 130
            self.esk = -90
            self.el_ichan2 = -60			#-60.06
            self.cao_ccanl = 2
            
            sec.insert("ccanl")
            for seg in sec:
                seg.catau_ccanl = 10		# Time constant for decay of intracellular Ca2+
                seg.caiinf_ccanl = 5.e-6		# Steady-state intracellular Ca2+ concentration

            
            sec.insert("borgka")
            for seg in sec:
                seg.gkabar_borgka = 0.00015		# A-type K+ conductance

            
            sec.insert("nca") # N-type Ca2+ conductance
            for seg in sec:
                seg.gncabar_nca = 0.0008   		# check to modify- original 0.004
           
            sec.insert("lca") 
            for seg in sec:
                seg.glcabar_lca = 0.005		# L-type Ca2+ conductance
           
            sec.insert("gskch") 
            for seg in sec:
                seg.gskbar_gskch = 0.000002		# Ca2+-dependent K (SK) conductance
           
            sec.insert("mykca") 
            for seg in sec:
                seg.gkbar_mykca = 0.0002			# Ca2+ and Voltage-dependent K+ (BK) conductance
		 					# make catau slower70e-3 	cao=2 cai=50.e-6





    def connect2target(self,target, delay = 1, weight=0.04): # { localobj nc #$o1 target point process, optional $o2 returned NetCon
        self.nc.append(h.NetCon(self.soma(0.5)._ref_v, target, sec=self.soma))
        self.nc[-1].threshold = -10 # mV
        self.nc[-1].delay = delay # ms
        self.nc[-1].weight[0] = weight # NetCon weight is a vector    
        return self.nc[-1]


    def addSynapses(self):
        self.pre_list = []
       
        # E0
        syn_ = h.MyExp2Syn(self.radM1(0.5))
        self.pre_list.append(syn_)    # AMPA        EC (not used)
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E1
        syn_ = h.MyExp2Syn(self.radM2(0.5))
        self.pre_list.append(syn_)    # AMPA        EC
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E2
        syn_ = h.MyExp2Syn(self.radM1(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E3
        syn_ = h.MyExp2Syn(self.radM2(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E4
        syn_ = h.MyExp2Syn(self.radT1(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E5
        syn_ = h.MyExp2Syn(self.radT2(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E6
        syn_ = h.MyExp2Syn(self.oriT1(0.5))
        self.pre_list.append(syn_)    # AMPA        PC
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E7
        syn_ = h.MyExp2Syn(self.oriT2(0.5))
        self.pre_list.append(syn_)    # AMPA        PC
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # I8
        syn_ = h.MyExp2Syn(self.soma(0.5))
        self.pre_list.append(syn_)    # GABA-A	Neighboring bistratified cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I9
        syn_ = h.MyExp2Syn(self.soma(0.5))
        self.pre_list.append(syn_)    # GABA-A	Basket cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I10
        syn_ = h.MyExp2Syn(self.oriT1(0.6))
        self.pre_list.append(syn_)    # GABA-A	Septum
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I11
        syn_ = h.MyExp2Syn(self.oriT2(0.6))
        self.pre_list.append(syn_)    # GABA-A	Septum
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I12
        syn_ = h.MyExp2Syn(self.oriT1(0.6))
        self.pre_list.append(syn_)    # GABA-B	Septum
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75

        
        # I13
        syn_ = h.MyExp2Syn(self.oriT2(0.6))
        self.pre_list.append(syn_)    # GABA-B	Septum
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75