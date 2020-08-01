# Dummy cell containing a BurstStim object
# BPG 10-12-08

class BurstCell():    
    def __init__(self):
        self.is_art =1
        self.ncstim = []
        self.stim = h.BurstStim2()
       	self.stim.number = 10000
       	self.stim.start = 0
       	self.stim.interval = 10
       	self.stim.noise = 0
       	self.stim.burstint = 100	# interburst interval (ms)
       	self.stim.burstlen = 100	# burst length (ms)

    def connect2target(self,target, delay = 1, weight=0.04): # { localobj nc #$o1 target point process, optional $o2 returned NetCon
        self.ncstim.append(h.NetCon(self.stim, target))
        self.ncstim[-1].delay = delay # ms
        self.ncstim[-1].weight[0] = weight # NetCon weight is a vector    
        return self.ncstim[-1]


