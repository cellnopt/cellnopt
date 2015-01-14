from cno.core.base import CNOBase

import pandas as pd
import numpy as np
import pylab
from cno.misc.profiler import do_profile
from cno.io.reactions import Reaction
import time


class Steady(CNOBase):
    """Naive implementation of Steady state to help in 
    designing the API"""
    
    def __init__(self, pknmodel, data, verbose=True):
        super(Steady, self).__init__(pknmodel, data, verbose)

        self.model = self.pknmodel.copy()
        # to speed up code later on
        self.model.buffer_reactions = self.model.reactions
        self.time = self.data.times[1]

        # just a reference to the conditions
        self.inhibitors = self.data.inhibitors
        self.stimuli = self.data.stimuli

        self.inhibitors_names = self.data.inhibitors.columns
        self.stimuli_names = self.data.stimuli.columns

        self.N = self.data.df.query('time==0').shape[0]

        self.results = self.data.df.copy()
        self.results = self.results.query("time==@self.time")

        # ignore data of the edge [0:2]
        self.toflip = [x[0:2] for x in self.model.edges(data=True) if x[2]['link'] == '-']

        self.init(self.time)

        self.measures = {}
        self.measures[0] = self.data.df.query("time==0").reset_index(drop=True).values
        self.measures[self.time] = self.data.df.query("time==@self.time").reset_index(drop=True).values

        self.simulated = {}

    def init(self, time):
        assert time in self.data.times
        self.values = {}
        for node in self.model.nodes():
            self.values[node] = np.array([np.nan for x in range(0,self.N)])

        for this in self.stimuli_names:
            self.values[this] = self.stimuli[this].values.copy()

        for this in self.inhibitors_names:
            self.values[this] = 1. - self.inhibitors[this].values.copy()

        self.and_gates = [x for x in self.model.nodes() if "^" in x]

        self.predecessors = {}
        for node in self.model.nodes():
            self.predecessors[node] = self.model.predecessors(node)

        self.successors = {}
        for node in self.model.nodes():
            self.successors[node] = self.model.successors(node)

    def preprocessing(self, expansion=True, compression=True, cutnonc=True):
        self.model.midas = self.data
        self.model.preprocessing(expansion=expansion, compression=compression, 
                cutnonc=cutnonc)
        self.init(self.time)
        self.model.buffer_reactions = self.model.reactions

    def reactions_to_predecessors(self, reactions):
    
        predecessors = dict(((k, []) for k in self.predecessors.keys()))
        for r in reactions:
            reac = Reaction(r)
            if "^" in reac.lhs:
                predecessors[r].extend(reac.lhs_species)
            else:
                predecessors[reac.rhs].append(reac.lhs.replace("!",""))
        return predecessors
            

    #@do_profile()
    def simulate(self, tick=1, debug=False, reactions=None):
        # pandas is very convenient but slower than numpy
        # The dataFrame instanciation is costly as well.
        # For small models, it has a non-negligeable cost.

        # inhibitors will be changed in not ON
        self.tochange = [x for x in self.model.nodes() if x not in self.stimuli_names 
                    and x not in self.and_gates ]

        # what about a species that is both inhibited and measured
        testVal = 1e-3

        values = self.values.copy()

        self.debug_values = []
        #self.X0 = pd.DataFrame(self.values)
        #self.debug_values.append(self.X0.copy())
        self.residuals = []
        self.penalties = []

        self.count = 0
        self.nSp = len(values.keys())
        residual = 1

        frac = 1.2
        # #FIXME +1 is to have same resrults as in CellnOptR
        # It means that if due to the cycles, you may not end up with same results.
        # this happends if you have cyvles with inhbititions
        # and an odd number of edges. 
        if reactions is None:
            # takes some time:
            # reactions = self.model.reactions
            reactions = self.model.buffer_reactions
            predecessors = self.predecessors.copy()
        else:
            predecessors = self.reactions_to_predecessors(reactions)
       
        values = self.values.copy()
        for inh in self.inhibitors:
            if len(predecessors[inh]) == 0:
                values[inh] = np.array([np.nan for x in range(0,self.N)])

 
        while (self.count < self.nSp * frac +1) and residual > testVal: 
            self.previous = values.copy()
            #self.X0 = pd.DataFrame(self.values)
            #self.X0 = self.values.copy()
            # compute AND gates first. why
            for node in self.and_gates:
                # replace na by large number so that min is unchanged
                if len(predecessors[node]) != 0:
                    values[node] = np.nanmin(np.array([values[x].copy() for x in 
                        predecessors[node]]), axis=0)
                else:
                    values[node] = self.previous[node]

            for node in self.tochange:
                # easy one, just the value of predecessors
                #if len(self.predecessors[node]) == 1:
                #    self.values[node] = self.values[self.predecessors[node][0]].copy()
                if len(predecessors[node]) == 0:
                    pass # nothing to change
                else:
                    values[node] = np.nanmax(
                            np.array([values[x] if (x,node) not in self.toflip 
                        else 1. - values[x] for x in 
                        predecessors[node]]),
                            axis=0)
                # take inhibitors into account
                if node in self.inhibitors_names:
                    values[node] *= 1 - self.inhibitors[node].values
            # 30 % of the time is here 
            # here NAs are set automatically to zero because of the int16 cast
            # but it helps speeding up a nit the code by removig needs to take care
            # of NAs. if we use sumna, na are ignored even when 1 is compared to NA            
            self.m1 = np.array([self.previous[k] for k in self.previous.keys() ], dtype=np.int16)
            self.m2 = np.array([values[k] for k in self.previous.keys() ], dtype=np.int16)
            residual = np.nansum(np.square(self.m1 - self.m2))

            self.debug_values.append(self.previous.copy())
           
            self.residuals.append(residual)
            self.count += 1

        # add the latest values simulated in the while loop
        self.debug_values.append(values.copy())

        # Need to set undefined values to NAs

        self.simulated[self.time] = np.array([values[k] 
            for k in self.data.df.columns ], dtype=float).transpose()
        self.prev = {}
        self.prev[self.time] = np.array([self.previous[k] 
            for k in self.data.df.columns ], dtype=float).transpose()

        mask = self.prev[self.time] != self.simulated[self.time]
        self.simulated[self.time][mask] = np.nan

        # set the non-resolved bits to NA
        # TODO
        #newInput[which(abs(outputPrev-newInput) > testVal)] <- NA

        # loops are handle diffenty

    #@do_profile()
    def score(self):
        # We need also to include NAFac, number of reactions in the model
        # for the sizeFac

        # on ExtLiverPCB
        # computeScoreT1(cnolist, pknmodel, rep(1,116))
        # 0.29199 from cellnoptr
        # found 0.2945
        # liverDREAM
        #computeScoreT1(cnolist, pknmodel, rep(1,58))
        #[1] 0.2574948
        # found 0.27
        # 0.2902250932121212

        # time 1 only is taken into account
        diff = np.square(self.measures[self.time] - self.simulated[self.time])
        #debugging
        self.diff = diff
        N = diff.shape[0] * diff.shape[1]
        Nna = np.isnan(diff).sum()
        N-= Nna
        #print(N)
        S = np.nansum(self.diff) / float(N)
        return S

    def plot(self):
        self.model.plot()

    def test(self, N=100):
        # N = 100, all bits on
        # CellNOptR on LiverDREAM  0.65 seconds. 0.9 in cno
        # CellNOptR on LiverDREAM preprocessed) on:0.75 seconds. 1.42 in cno
        # 0.2574948 in CellNOptR
        # 0.27

        # CellNOptR on ToyMMB      : 0.22        ; 0.22s in cno
        # 0.09467
        # process and  "EGF=Raf"        "EGF+TNFa=PI3K"  "Erk+TNFa=Hsp27" off
        # then MSE is 0.10838/2

        # CellNOptR on ExtLiverPCB : 1.4 seconds ; 1.7s in cno
        # 0.29199

        # cost of pruning models ?
        """
        library(CellNOptR)
        cnolist = ...
        pknmodel = ...
        test <- function(cnolist, model, N=100){
        t1 = proc.time()
         for (i in rep(1,N)){ 
            mse = computeScoreT1(cnolist, model, rep(1,length(model$reacID)))
         }
         t2 = proc.time()
         print(t2-t1)
         return(mse)
         }
         test(cnolist, model)
        """
        t1 = time.time()
        for i in range(0,N):
            self.init(self.time)
            self.simulate()
            self.score()
        t2 = time.time()
        print(str(t2-t1) + " seconds")        

    def plotsim(self, fontsize=16):
        cm = pylab.get_cmap('gray')
        pylab.clf()
        data = pd.DataFrame(self.debug_values[-1]).fillna(0.5)
        pylab.pcolor(data, cmap=cm, vmin=0, vmax=1,
                shading='faceted')
        pylab.colorbar()
        ax1 = pylab.gca()
        ax1.set_xticks([])
        Ndata = len(data.columns)
        ax1.set_xlim(0, Ndata)
        ax = pylab.twiny()
        ax.set_xticks(pylab.linspace(0.5, Ndata+0.5, Ndata ))
        ax.set_xticklabels(data.columns, fontsize=fontsize, rotation=90)
        times = list(data.index)
        Ntimes = len(times)
        ax1.set_yticks([x+0.5 for x in times])
        ax1.set_yticklabels(times[::-1],
                    fontsize=fontsize)
        pylab.sca(ax1)
        pylab.title("Steady state for all experiments(x-axis)\n\n\n\n")
        pylab.tight_layout()

    def plot_errors(self, columns=None):
        # What do we use here: self.values
        print("Use only time 0..")

        if columns is None:
            columns = self.data.df.columns
        X1 = pd.DataFrame(self.debug_values[-1])[columns].copy()
        N = X1.shape[0]

        X1['time'] = [self.time] * N
        X1['cell'] = [self.data.cellLine] * N
        X1['experiment'] = self.data.experiments.index
        X1.set_index(['cell', 'experiment', 'time'], inplace=True)
        self.data.sim.ix[X1.index] = X1
        self.data.plot(mode='mse')
        print("MSE= %s(caspo/cno with only 1 time)" % self.score())
        print("MSE= %s(cellnoptr with only 1 time)" % str(self.score()/2.))



#@do_profile()
def test():
    from cno import cnodata
    s = Steady(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
    s.test()
#test()







"""

$mse
[,1]  [,2]   [,3] [,4]   [,5]    [,6] [,7]
[1,] 0.00405 0.500 0.3698 0.02 0.0072 0.00000 0.00
[2,] 0.01620 0.045 0.0050 0.00 0.0000 0.28125 0.18
[3,] 0.00405 0.045 0.0050 0.02 0.0072 0.28125 0.18
[4,] 0.00405 0.000 0.3698 0.00 0.0000 0.00000 0.00
[5,] 0.01620 0.045 0.0050 0.00 0.0000 0.28125 0.18
[6,] 0.00405 0.045 0.0050 0.00 0.0000 0.28125 0.18
[7,] 0.00000 0.500 0.0000 0.02 0.0072 0.00000 0.00
[8,] 0.00000 0.045 0.0050 0.50 0.5000 0.28125 0.18
[9,] 0.00000 0.045 0.0050 0.02 0.0072 0.28125 0.18

$simResults[[1]]$t0
[,1] [,2] [,3] [,4] [,5] [,6] [,7]
[1,]    0    0    0    0    0    0    0
[2,]    0    0    0    0    0    0    0
[3,]    0    0    0    0    0    0    0
[4,]    0    0    0    0    0    0    0
[5,]    0    0    0    0    0    0    0
[6,]    0    0    0    0    0    0    0
[7,]    0    0    0    0    0    0    0
[8,]    0    0    0    0    0    0    0
[9,]    0    0    0    0    0    0    0

s.score()
pd.DataFrame(s.diff/2.)[[0,2,4,1,6,3,5]]

     Akt Hsp27 NFkB Erk p90RSK Jnk cJun
[,1] [,2] [,3] [,4] [,5] [,6] [,7]
[1,]    1    1    0    1    1    0    0
[2,]    1    1    1    0    0    1    1
[3,]    1    1    1    1    1    1    1
[4,]    1    0    0    0    0    0    0
[5,]    1    1    1    0    0    1    1
[6,]    1    1    1    0    0    1    1
[7,]    0    1    0    1    1    0    0
[8,]    0    1    1    1    1    1    1
[9,]    0    1    1    1    1    1    1

    EGF TNFa Raf PI3K
[1,]   1    0   0    0
[2,]   0    1   0    0
[3,]   1    1   0    0
[4,]   1    0   1    0
[5,]   0    1   1    0
[6,]   1    1   1    0
[7,]   1    0   0    1
[8,]   0    1   0    1
[9,]   1    1   0    1

"""


"""

$mse
[,1]       [,2]      [,3]      [,4]      [,5]       [,6]
[1,] 0.02607682         NA 0.3539345 0.4100143        NA 0.03221782
[2,] 0.09200486 0.03723420 0.2499203 0.6794323 0.3669156 0.02610956
[3,] 0.02845594         NA 0.1943445 0.7901995 0.3425479 0.28019457
[4,] 0.07943784 0.04197267 0.2153080 0.5394515 0.3405881 0.26677940
[5,] 0.06059407 0.03067175 0.2714396 0.4428013 0.3980180 0.03146179
[6,] 0.02276370         NA 0.2004797 0.3669970 0.3474665 0.28991532
[7,] 0.06374754 0.02499234 0.2358515 0.4051598 0.4542814 0.20409221
[8,] 0.01885474 0.01503284 0.2139463 0.6602201 0.3880110 0.03740429
[9,] 0.01601736 0.02506731 0.2179543 0.7714214 0.4299402 0.24682266
[10,] 0.02014269 0.01944989 0.2992465 0.6237731 0.3529949 0.25041879

$simResults
$simResults[[1]]
$simResults[[1]]$t0
[,1] [,2] [,3] [,4] [,5] [,6]
[1,]    0    0    0    0    0    0
[2,]    0    0    0    0    0    0
[3,]    0    0    0    0    0    0
[4,]    0    0    0    0    0    0
[5,]    0    0    0    0    0    0
[6,]    0    0    0    0    0    0
[7,]    0    0    0    0    0    0
[8,]    0    0    0    0    0    0
[9,]    0    0    0    0    0    0
[10,]    0    0    0    0    0    0

$simResults[[1]]$t1
[,1] [,2] [,3] [,4] [,5] [,6]
[1,]    0   NA    1    1   NA    0
[2,]    1    1    1    0    1    0
[3,]    0   NA    1    0    1    0
[4,]    1    1    1    0    1    0
[5,]    1    1    1    1    1    0
[6,]    0   NA    1    1    1    0
[7,]    1    1    1    1    1    0
[8,]    0    0    1    0    1    0
[9,]    0    0    1    0    1    0
[10,]    0    0    1    0    1    0


imResults[[1]]$t1
      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
 [1,]    1    0    0    0    0    0    0
 [2,]    1    1    1    0    0    1    1
 [3,]    1    1    1    0    0    1    1
 [4,]    1    0    0    0    0    0    0
 [5,]    1    1    1    0    0    1    1
 [6,]    1    1    1    0    0    1    1
 [7,]    0    1    0    1    1    0    0
 [8,]    0    1    1    1    1    1    1
 [9,]    0    1    1    1    1    1    1






"""
