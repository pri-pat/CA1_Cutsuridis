# Data from Saraga et al. (2003) paper
# changed morphology and some channel densities (BPG 12-1-09)
#   OLM_Cell
# basic_shape,
# pre_list, connect2target
# soma, dend1, dend2, axon
# all
from neuron import h

class PyramidalCell():
    """ Pyramidal Cell definition """
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
        return "Pyramidal Cell {}".format(self.gid)

    def create_sections(self):
        self.soma = h.Section(name='soma', cell=self)
        self.radTprox = h.Section(name='radTprox', cell=self)
        self.radTmed = h.Section(name='radTmed', cell=self)
        self.radTdist = h.Section(name='radTdist', cell=self)
        self.axon = h.Section(name='axon', cell=self)
        self.lm_thick2 = h.Section(name='lm_thick2', cell=self)
        self.lm_thick1 = h.Section(name='lm_thick1', cell=self)
        self.oriprox1 = h.Section(name='oriprox1', cell=self)
        self.oridist1 = h.Section(name='oridist1', cell=self)
        self.oriprox2 = h.Section(name='oriprox2', cell=self)
        self.oridist2 = h.Section(name='oridist2', cell=self)
        self.lm_medium2 = h.Section(name='lm_medium2', cell=self)
        self.lm_thin2 = h.Section(name='lm_thin2', cell=self)
        self.lm_medium1 = h.Section(name='lm_medium1', cell=self)
        self.lm_thin1 = h.Section(name='lm_thin1', cell=self)

    def build_topology(self):
        self.radTprox.connect(self.soma(1))
        self.radTmed.connect(self.radTprox(1))
        self.radTdist.connect(self.radTmed(1))
        self.axon.connect(self.soma(1))
        self.lm_thick2.connect(self.radTdist(1))
        self.lm_medium2.connect(self.lm_thick2(1))
        self.lm_thin2.connect(self.lm_medium2(1))

        self.lm_thick1.connect(self.radTdist(1))
        self.lm_medium1.connect(self.lm_thick1(1))
        self.lm_thin1.connect(self.lm_medium1(1))

        self.oriprox1.connect(self.soma(0))
        self.oridist1.connect(self.oriprox1(1))

        self.oriprox2.connect(self.soma(1))
        self.oridist2.connect(self.oriprox2(1))
      

    def define_geometry(self):
        for sec in self.all:
            sec.pt3dclear()

        h.pt3dadd(0, 0, 0, 10, sec=self.soma)
        h.pt3dadd(15, 0, 0, 10, sec=self.soma)
        
        h.pt3dadd(15, 0, 0, 1, sec=self.radTprox)
        h.pt3dadd(15, 30, 0, 1, sec=self.radTprox)
       
        h.pt3dadd(15, 30, 0, 1, sec=self.radTmed)
        h.pt3dadd(15, 60, 0, 1, sec=self.radTmed)
        
        h.pt3dadd(15, 60, 0, 1, sec=self.radTdist)
        h.pt3dadd(15, 90, 0, 1, sec=self.radTdist)

        h.pt3dadd(15, 90, 0, 1, sec=self.lm_thick2)
        h.pt3dadd(45, 105, 0, 1, sec=self.lm_thick2)

        h.pt3dadd(45, 105, 0, 1, sec=self.lm_medium2)
        h.pt3dadd(75, 120, 0, 1, sec=self.lm_medium2)

        h.pt3dadd(75, 120, 0, 1, sec=self.lm_thin2)
        h.pt3dadd(105, 135, 0, 1, sec=self.lm_thin2)

        h.pt3dadd(15, 90, 0, 1, sec=self.lm_thick1)
        h.pt3dadd(-14, 105, 0, 1, sec=self.lm_thick1)

        h.pt3dadd(-14, 105, 0, 1, sec=self.lm_medium1)
        h.pt3dadd(-44, 120, 0, 1, sec=self.lm_medium1)

        h.pt3dadd(-44, 120, 0, 1, sec=self.lm_thin1)
        h.pt3dadd(-89, 135, 0, 1, sec=self.lm_thin1)
         
        h.pt3dadd(15, 0, 0, 1, sec=self.axon)
        h.pt3dadd(15, -149, 0, 1, sec=self.axon)

        h.pt3dadd(0, 0, 0, 1, sec=self.oriprox1)
        h.pt3dadd(-44, -29, 0, 1, sec=self.oriprox1)

        h.pt3dadd(-44, -29, 0, 1, sec=self.oridist1)
        h.pt3dadd(-74, -59, 0, 1, sec=self.oridist1)

        h.pt3dadd(15, 0, 0, 1, sec=self.oriprox2)
        h.pt3dadd(60, -29, 0, 1, sec=self.oriprox2)

        h.pt3dadd(60, -29, 0, 1, sec=self.oridist2)
        h.pt3dadd(105, -59, 0, 1, sec=self.oridist2)

        
        self.soma.L = 20
        self.soma.diam = 10
        self.radTprox.L = 100
        self.radTprox.diam = 4
        self.radTmed.L = 100
        self.radTmed.diam = 3
        self.radTdist.L = 200
        self.radTdist.diam = 2

        self.lm_medium2.L = 100
        self.lm_medium2.diam = 1.5
        self.lm_thin2.L = 100
        self.lm_thin2.diam = 1
        self.lm_thick2.L = 100
        self.lm_thick2.diam = 2
        
        
        self.lm_thick1.L = 200
        self.lm_thick1.diam = 2
        self.lm_medium1.L = 100
        self.lm_medium1.diam = 1.5
        self.lm_thin1.L = 100
        self.lm_thin1.diam = 1

        self.oriprox1.L = 100
        self.oriprox1.diam = 2
        self.oridist1.L = 200
        self.oridist1.diam = 1.5

        self.oriprox2.L = 100
        self.oriprox2.diam = 2
        self.oridist2.L = 200
        self.oridist2.diam = 1.5


        self.axon.L = 100
        self.axon.diam = 1
 
        for sec in self.all:
            h("lf = lambda_f(100)")
            sec.nseg = int((sec.L/(0.1*h.lf)+.9)/2)*2 + 1 

    def build_subsets(self):
        self.all = h.SectionList()
        self.all.wholetree(sec=self.soma)


    def define_biophysics(self):

        gka_soma = 0.0075
        gh_soma = 0.00005
        Rm = 20000    # 28000 Ohm.cm^2 (Migliore value)
        
        self.soma.insert("hha2") # HH mechanism with low threshold for Na spikes (-57 mV)
        self.soma.insert("pas") # leak conductance
        self.soma.insert("h") # h current according to Migliore et al. 2004 
        self.soma.insert("kap") # proximal A current
        self.soma.insert("km") # m-type potassium current
        self.soma.insert("cal") # HVA Ca++-L type current
        self.soma.insert("cat") # LVA Ca++-T type current
        self.soma.insert("somacar") # HVAm Ca++-R type current
        self.soma.insert("kca") # K(Ca) sAHP potassium type current
        self.soma.insert("mykca") # medium AHP K++ current (BPG)
        self.soma.insert("cad") # calcium pump/buffering mechanism
        for seg in self.soma:
            seg.gnabar_hha2 = 0.007
            seg.gkbar_hha2  = 0.007/5
            seg.gl_hha2     = 0
            seg.el_hha2     = -70
            seg.g_pas =  1/Rm

            seg.ghdbar_h = gh_soma
            seg.vhalfl_h = -73
            
            seg.gkabar_kap = gka_soma        # 0.0075
            seg.gbar_km    = 0.06
            seg.gcalbar_cal = 0.0014/2
            seg.gcatbar_cat = 0.0001/2
            seg.gcabar_somacar = 0.0003
            seg.gbar_kca = 5*0.0001
            seg.gkbar_mykca = 0.09075

          
  
  # //        insert hNa            // h current according to Poirazi 2003
  # //        gbar_h  = 0.000043        // anything above 0.000043 gives hyperpolarizing oscillations
  # //        gbar_h  = 1.872e-5        
  # //        K_h     = 8.8
  # //        vhalf_h = -82
  
        self.radTprox.insert("h") # h current according to Migliore et al. 2004 
        self.radTprox.insert("car") # h current according to Migliore et al. 2004 
        self.radTprox.insert("calH") # h current according to Migliore et al. 2004 
        self.radTprox.insert("cat") # h current according to Migliore et al. 2004 
        self.radTprox.insert("cad") # calcium pump/buffering mechanism
        self.radTprox.insert("kca") # slow AHP K+ current
        self.radTprox.insert("mykca") # medium AHP K++ current (BPG)
        self.radTprox.insert("km") # m-type K current
        self.radTprox.insert("kap") # Inserting A-current
        self.radTprox.insert("kad") # Inserting A-current
        self.radTprox.insert("hha_old") # HH mechanism with high threshold for Na spikes (-50 mV)
        self.radTprox.insert("pas") # passive
        for seg in self.radTprox:
            seg.ghdbar_h = 2*gh_soma # 0.000005    
            seg.vhalfl_h = -81
            seg.gcabar_car = 0.1*0.0003
            seg.gcalbar_calH = 0.1*0.00031635    # varies from .1*0.00031635 to 4.6*0.00031635 as distance increases
            seg.gcatbar_cat = 0.0001
            seg.gbar_kca = 5*0.0001        # varies depending on distance from 0.5*0.0001 to 5*0.0001
            seg.gkbar_mykca = 2*0.0165
            seg.gbar_km = 0.06        # varies with distance (see Poirazzi et al. 2003 cell-setup.hoc file)
            seg.gkabar_kap = 2*gka_soma        # 0.0075
            seg.gkabar_kad = 0
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old     = -70
        
# //        insert hNa            // h current according to Poirazi 2003
# //              gbar_h  = 0.000043        // anything above 0.000043 gives hyperpolarizing oscillations
# //             gbar_h  = 1.872e-5        
# //                K_h     = 8.8
# //                vhalf_h = -82

        self.radTmed.insert("h") # h current according to Migliore et al. 2004
        self.radTmed.insert("car") # HVAm Ca++-R type current
        self.radTmed.insert("calH") # HVA L-type Ca2+ channel used in distal dendrites to account for distally restricted initiation of Ca2+ spikes
        self.radTmed.insert("cat") # HVA T-type Ca2+ channel 
        self.radTmed.insert("cad") # calcium pump/buffering mechanism
        self.radTmed.insert("kca") # slow AHP K+ current
        self.radTmed.insert("mykca") # medium AHP K++ current (BPG)
        self.radTmed.insert("km") # m-type K current
        self.radTmed.insert("kap") # Inserting A-current
        self.radTmed.insert("kad") 
        self.radTmed.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.radTmed.insert("pas") # leak conductance
        
        for seg in self.radTmed:
            seg.ghdbar_h = 4*gh_soma            # 0.000005                    
            seg.vhalfl_h = -81
            seg.gcabar_car = 0.1*0.0003
            seg.gcalbar_calH = 10*0.00031635    # 4.6*0.00031635 varies from .1*0.00031635 to 4.6*0.00031635 as distance increases
            seg.gcatbar_cat = 0.0001        # 0.0001
            seg.gbar_kca = 5*0.0001        # varies depending on distance from 0.5*0.0001 to 5*0.0001
            seg.gkbar_mykca = 2*0.0165
            seg.gbar_km = 0.06            # varies with distance (see Poirazzi et al. 2003 cell-setup.hoc file)
            seg.gkabar_kap = 0
            seg.gkabar_kad = 4*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old     = -70
            
# //        insert hNa            // h current according to Poirazi 2003
# //              gbar_h  = 0.000043        
# //             gbar_h  = 1.872e-5        
# //                K_h     = 8.8
# //                vhalf_h = -82


        self.radTdist.insert("h") # h current according to Migliore et al. 2004
        self.radTdist.insert("car") # HVAm Ca++-R type current
        self.radTdist.insert("calH") # HVA L-type Ca2+ channel used in distal dendrites to account for distally restricted initiation of Ca2+ spikes
        self.radTdist.insert("cat") # HVA T-type Ca2+ channel 
        self.radTdist.insert("cad") # calcium pump/buffering mechanism
        self.radTdist.insert("kca") # slow AHP K+ current
        self.radTdist.insert("mykca") # medium AHP K++ current (BPG)
        self.radTdist.insert("km") # m-type K current
        self.radTdist.insert("kap") # Inserting A-current
        self.radTdist.insert("kad") 
        self.radTdist.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.radTdist.insert("pas") # leak conductance
        
        for seg in self.radTdist:
            seg.ghdbar_h = 7*gh_soma            # 0.000005                    
            seg.vhalfl_h = -81
            seg.gcabar_car = 0.1*0.0003
            seg.gcalbar_calH = 10*0.00031635    # 4.6*0.00031635 varies from .1*0.00031635 to 4.6*0.00031635 as distance increases
            seg.gcatbar_cat = 0.0001        # 0.0001
            seg.gbar_kca = 0.5*0.0001        # varies depending on distance from 0.5*0.0001 to 5*0.0001
            seg.gkbar_mykca = 0.25*0.0165
            seg.gbar_km = 0.06            # varies with distance (see Poirazzi et al. 2003 cell-setup.hoc file)
            seg.gkabar_kap = 0
            seg.gkabar_kad = 6*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old = -70



        self.lm_thick2.insert("kad") 
        self.lm_thick2.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.lm_thick2.insert("pas") # leak conductance
        
        for seg in self.lm_thick2:
            seg.gkabar_kad = 6.5*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old = -70
            seg.g_pas = 1/200000


        self.lm_medium2.insert("kad") 
        self.lm_medium2.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.lm_medium2.insert("pas") # leak conductance
        
        for seg in self.lm_medium2:
            seg.gkabar_kad = 6.5*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old = -70
            seg.g_pas = 1/200000


        self.lm_thin2.insert("kad") 
        self.lm_thin2.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.lm_thin2.insert("pas") # leak conductance
        
        for seg in self.lm_thin2:
            seg.gkabar_kad = 6.5*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old = -70
            seg.g_pas = 1/200000


        self.lm_thick1.insert("kad") 
        self.lm_thick1.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.lm_thick1.insert("pas") # leak conductance
        
        for seg in self.lm_thick1:
            seg.gkabar_kad = 6.5*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old = -70
            seg.g_pas = 1/200000


        self.lm_medium1.insert("kad") 
        self.lm_medium1.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.lm_medium1.insert("pas") # leak conductance
        
        for seg in self.lm_medium1:
            seg.gkabar_kad = 6.5*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old = -70
            seg.g_pas = 1/200000


        self.lm_thin1.insert("kad") 
        self.lm_thin1.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.lm_thin1.insert("pas") # leak conductance
        
        for seg in self.lm_thin1:
            seg.gkabar_kad = 6.5*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old = -70
            seg.g_pas = 1/200000



        self.oriprox1.insert("h") # h current according to Migliore et al. 2004
        self.oriprox1.insert("car") # HVAm Ca++-R type current
        self.oriprox1.insert("calH") # HVA L-type Ca2+ channel used in distal dendrites to account for distally restricted initiation of Ca2+ spikes
        self.oriprox1.insert("cat") # HVA T-type Ca2+ channel 
        self.oriprox1.insert("cad") # calcium pump/buffering mechanism
        self.oriprox1.insert("kca") # slow AHP K+ current
        self.oriprox1.insert("mykca") # medium AHP K++ current (BPG)
        self.oriprox1.insert("km") # m-type K current
        self.oriprox1.insert("kap") # Inserting A-current
        self.oriprox1.insert("kad") 
        self.oriprox1.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.oriprox1.insert("pas") # leak conductance
        
        for seg in self.oriprox1:
            seg.ghdbar_h = gh_soma            # 0.000005                    
            seg.vhalfl_h = -81
            seg.gcabar_car = 0.1*0.0003
            seg.gcalbar_calH = 0.1*0.00031635 # 4.6*0.00031635 varies from .1*0.00031635 to 4.6*0.00031635 as distance increases
            seg.gcatbar_cat = 0.0001        # 0.0001
            seg.gbar_kca = 5*0.0001        # varies depending on distance from 0.5*0.0001 to 5*0.0001
            seg.gkbar_mykca = 2*0.0165
            seg.gbar_km = 0.06            # varies with distance (see Poirazzi et al. 2003 cell-setup.hoc file)
            seg.gkabar_kap = gka_soma
            seg.gkabar_kad = 0
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old     = -70
            

        self.oriprox2.insert("h") # h current according to Migliore et al. 2004
        self.oriprox2.insert("car") # HVAm Ca++-R type current
        self.oriprox2.insert("calH") # HVA L-type Ca2+ channel used in distal dendrites to account for distally restricted initiation of Ca2+ spikes
        self.oriprox2.insert("cat") # HVA T-type Ca2+ channel 
        self.oriprox2.insert("cad") # calcium pump/buffering mechanism
        self.oriprox2.insert("kca") # slow AHP K+ current
        self.oriprox2.insert("mykca") # medium AHP K++ current (BPG)
        self.oriprox2.insert("km") # m-type K current
        self.oriprox2.insert("kap") # Inserting A-current
        self.oriprox2.insert("kad") 
        self.oriprox2.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.oriprox2.insert("pas") # leak conductance
        
        for seg in self.oriprox2:
            seg.ghdbar_h = gh_soma            # 0.000005                    
            seg.vhalfl_h = -81
            seg.gcabar_car = 0.1*0.0003
            seg.gcalbar_calH = 0.1*0.00031635 # 4.6*0.00031635 varies from .1*0.00031635 to 4.6*0.00031635 as distance increases
            seg.gcatbar_cat = 0.0001        # 0.0001
            seg.gbar_kca = 5*0.0001        # varies depending on distance from 0.5*0.0001 to 5*0.0001
            seg.gkbar_mykca = 2*0.0165
            seg.gbar_km = 0.06            # varies with distance (see Poirazzi et al. 2003 cell-setup.hoc file)
            seg.gkabar_kap = 0.0075
            seg.gkabar_kad = 0
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old     = -70


        self.oridist1.insert("h") # h current according to Migliore et al. 2004
        self.oridist1.insert("car") # HVAm Ca++-R type current
        self.oridist1.insert("calH") # HVA L-type Ca2+ channel used in distal dendrites to account for distally restricted initiation of Ca2+ spikes
        self.oridist1.insert("cat") # HVA T-type Ca2+ channel 
        self.oridist1.insert("cad") # calcium pump/buffering mechanism
        self.oridist1.insert("kca") # slow AHP K+ current
        self.oridist1.insert("mykca") # medium AHP K++ current (BPG)
        self.oridist1.insert("km") # m-type K current
        self.oridist1.insert("kap") # Inserting A-current
        self.oridist1.insert("kad") 
        self.oridist1.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.oridist1.insert("pas") # leak conductance
        
        for seg in self.oridist1:
            seg.ghdbar_h = 2*gh_soma            # 0.000005                    
            seg.vhalfl_h = -81
            seg.gcabar_car = 0.1*0.0003
            seg.gcalbar_calH = 0.1*0.00031635 # 4.6*0.00031635 varies from .1*0.00031635 to 4.6*0.00031635 as distance increases
            seg.gcatbar_cat = 0.0001        # 0.0001
            seg.gbar_kca = 5*0.0001        # varies depending on distance from 0.5*0.0001 to 5*0.0001
            seg.gkbar_mykca = 2*0.0165
            seg.gbar_km = 0.06            # varies with distance (see Poirazzi et al. 2003 cell-setup.hoc file)
            seg.gkabar_kap = gka_soma
            seg.gkabar_kad = 0
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old     = -70

        self.oridist2.insert("h") # h current according to Migliore et al. 2004
        self.oridist2.insert("car") # HVAm Ca++-R type current
        self.oridist2.insert("calH") # HVA L-type Ca2+ channel used in distal dendrites to account for distally restricted initiation of Ca2+ spikes
        self.oridist2.insert("cat") # HVA T-type Ca2+ channel 
        self.oridist2.insert("cad") # calcium pump/buffering mechanism
        self.oridist2.insert("kca") # slow AHP K+ current
        self.oridist2.insert("mykca") # medium AHP K++ current (BPG)
        self.oridist2.insert("km") # m-type K current
        self.oridist2.insert("kap") # Inserting A-current
        self.oridist2.insert("kad") 
        self.oridist2.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.oridist2.insert("pas") # leak conductance
        
        for seg in self.oridist2:
            seg.ghdbar_h = 2*gh_soma            # 0.000005                    
            seg.vhalfl_h = -81
            seg.gcabar_car = 0.1*0.0003
            seg.gcalbar_calH = 0.1*0.00031635 # 4.6*0.00031635 varies from .1*0.00031635 to 4.6*0.00031635 as distance increases
            seg.gcatbar_cat = 0.0001        # 0.0001
            seg.gbar_kca = 5*0.0001        # varies depending on distance from 0.5*0.0001 to 5*0.0001
            seg.gkbar_mykca = 2*0.0165
            seg.gbar_km = 0.06            # varies with distance (see Poirazzi et al. 2003 cell-setup.hoc file)
            seg.gkabar_kap = 0.0075
            seg.gkabar_kad = 0
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old     = -70



        self.axon.insert("km") 
        self.axon.insert("hha2") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.axon.insert("pas") # leak conductance
        
        for seg in self.axon:
            seg.gbar_km     = 0.5*0.06
            seg.gnabar_hha2 = 0.1
            seg.gkbar_hha2  = 0.1/5
            seg.gl_hha2 = 0
            seg.el_hha2 = -70
            seg.g_pas = 1/Rm


        for sec in self.all:
            # self.cm = Not setting cm
            sec.Ra = 150            # 31.3 +/- 10.9
            sec.cm = 1
            sec.ena = 50
            sec.e_pas = -70
            sec.g_pas = 1/Rm # crucial parameter for backpropagating action potential spiking of PCs
            sec.ek = -80
            


    def connect2target(self,target, delay = 1, weight=0.04): # { localobj nc #$o1 target point process, optional $o2 returned NetCon
        self.nc.append(h.NetCon(self.soma(0.5)._ref_v, target, sec=self.soma))
        self.nc[-1].threshold = -10 # mV
        self.nc[-1].delay = delay # ms
        self.nc[-1].weight[0] = weight # NetCon weight is a vector    
        return self.nc[-1]


    def addSynapses(self):
        self.pre_list = []

        # E0
        syn_ = h.MyExp2Syn(self.lm_thick1(0.5))
        self.pre_list.append(syn_)    # AMPA        EC
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # E1
        syn_ = h.MyExp2Syn(self.lm_thick2(0.5))
        self.pre_list.append(syn_)    # AMPA        EC 
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # E2
        syn_ = h.MyExp2Syn(self.radTmed(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E3
        syn_ = h.NMDA(self.radTmed(0.5))
        self.pre_list.append(syn_)    # NMDA        CA3 Shaffer collateral
        syn_.tcon = 2.3    
        syn_.tcoff = 100
        syn_.gNMDAmax = 1    # use connection weight to determine max cond

        # E4
        syn_ = h.MyExp2Syn(self.radTprox(0.5))
        self.pre_list.append(syn_)    # AMPA        PC Recurrent collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # I5
        syn_ = h.MyExp2Syn(self.soma(0.5))
        self.pre_list.append(syn_)    # GABA-A    basket cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I6
        syn_ = h.MyExp2Syn(self.axon(0.1))
        self.pre_list.append(syn_)    # GABA-A    AA cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I7
        syn_ = h.MyExp2Syn(self.lm_thick1(0.5))
        self.pre_list.append(syn_)    # GABA-A    OLM cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I8
        syn_ = h.MyExp2Syn(self.lm_thick2(0.5))
        self.pre_list.append(syn_)    # GABA-A    OLM cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I9
        syn_ = h.MyExp2Syn(self.lm_thick1(0.5))
        self.pre_list.append(syn_)    # GABA-B    OLM cell
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75
        
        # I10
        syn_ = h.MyExp2Syn(self.lm_thick2(0.5))
        self.pre_list.append(syn_)    # GABA-B    OLM Cell
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75
        
        # I11
        syn_ = h.MyExp2Syn(self.radTmed(0.8))
        self.pre_list.append(syn_)    # GABA-A    Bistratified
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I12
        syn_ = h.MyExp2Syn(self.radTmed(0.7))
        self.pre_list.append(syn_)    # GABA-A    Bistratified
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75

        
        # I13
        syn_ = h.MyExp2Syn(self.radTmed(0.6))
        self.pre_list.append(syn_)    # GABA-A    Bistratified
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
                
        # I14
        syn_ = h.MyExp2Syn(self.radTmed(0.4))
        self.pre_list.append(syn_)    # GABA-A    Bistratified
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
                
        # I15
        syn_ = h.MyExp2Syn(self.radTmed(0.3))
        self.pre_list.append(syn_)    # GABA-A    Bistratified
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75       
                
        # I16
        syn_ = h.MyExp2Syn(self.radTmed(0.2))
        self.pre_list.append(syn_)    # GABA-A    Bistratified
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75

                
        # I17
        syn_ = h.MyExp2Syn(self.radTmed(0.8))
        self.pre_list.append(syn_)    # GABA-B    Bistratified
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75
        
                
        # I18
        syn_ = h.MyExp2Syn(self.radTmed(0.7))
        self.pre_list.append(syn_)    # GABA-B    Bistratified
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75
        
                
        # I19
        syn_ = h.MyExp2Syn(self.radTmed(0.6))
        self.pre_list.append(syn_)    # GABA-B    Bistratified
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75        
                
        # I20
        syn_ = h.MyExp2Syn(self.radTmed(0.4))
        self.pre_list.append(syn_)    # GABA-B    Bistratified
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75

                
        # I21
        syn_ = h.MyExp2Syn(self.radTmed(0.3))
        self.pre_list.append(syn_)    # GABA-B    Bistratified
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75       
                
        # I22
        syn_ = h.MyExp2Syn(self.radTmed(0.2))
        self.pre_list.append(syn_)    # GABA-B    Bistratified
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75

                
        # I23
        syn_ = h.STDPE2(self.radTmed(0.5))
        self.pre_list.append(syn_)    # AMPA modifiable	CA3 Schaffer collaterals
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0