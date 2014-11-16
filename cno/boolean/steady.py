from cno.core.base import CNOBase
import pandas as pd
import numpy as np
import pylab


class Steady(CNOBase):
    """Naive implementation of Steady state to help in 
    designing the API"""
    
    def __init__(self, model, data, verbose=True):
        super(Steady, self).__init__(model, data, verbose)

        self.model = self.pknmodel.copy()
        self.time = self.data.times[1]

        # just a reference to the conditions
        self.inhibitors = self.data.inhibitors
        self.stimuli = self.data.stimuli

        self.inhibitors_names = self.data.inhibitors.columns
        self.stimuli_names = self.data.stimuli.columns

        self.predecessors = {}
        for node in self.model.nodes():
            self.predecessors[node] = self.model.predecessors(node)

        self.successors = {}
        for node in self.model.nodes():
            self.successors[node] = self.model.successors(node)

        self.and_gates = [x for x in self.model.nodes() if "^" in x]
        self.results = self.data.df.copy()
        self.results = self.results.query("time==@self.time")

        # ignore data of the edge [0:2]
        self.toflip = [x[0:2] for x in self.model.edges(data=True) if x[2]['link'] == '-']

        self.init(self.time)

    def init(self, time):
        assert time in self.data.times
        # this should have the same size as experiments ?
        # FIXME is it correct
        N = self.data.df.query('time==@time').shape[0]

        self.values = {}
        for node in self.model.nodes():
            self.values[node] = np.array([np.nan for x in range(0,N)])

        for this in self.stimuli_names:
            self.values[this] = self.stimuli[this].values.copy()
        for this in self.inhibitors_names:
            self.values[this] = 1. - self.inhibitors[this].values.copy()

    def preprocessing(self):
        self.model.midas = self.data
        self.model.preprocessing()
        self.init(self.time)

    def simulate(self, tick, debug=False):
        #self.tochange = [x for x in self.model.nodes() if x not in self.stimuli_names 
        #            and x not in self.inhibitors_names and x not in self.and_gates ]
        self.tochange = [x for x in self.model.nodes() if x not in self.stimuli_names 
                    and x not in self.and_gates ]

        # what about a species that is both inhibited and measured
        count = 0
        testVal = 1e-3
        nSp = max(2, len(self.data.df.columns))
        residual = 1

        # FIXME: use numpy.matrix instead of dataframe

        self.debug_values = []
        self.X0 = pd.DataFrame(self.values)
        self.debug_values.append(self.X0.copy())
        self.residuals = []


        while (count < nSp * 1.2 ) and residual > testVal: 
            self.X0 = pd.DataFrame(self.values)
            #self.X0 = self.values.copy()
            # compute AND gates first. why
            for node in self.and_gates:
                self.values[node] = np.array([self.values[x].copy() for x in 
                        self.predecessors[node]]).min(axis=0)

            for node in self.tochange:
                # easy one, just the value of predecessors
                if len(self.predecessors[node]) == 1:
                    self.values[node] = self.values[self.predecessors[node][0]].copy()
                elif len(self.predecessors[node]) == 0:
                    pass # nothing to change
                else:
                    self.values[node] = np.array([self.values[x] if (x,node) not in self.toflip 
                        else 1-self.values[x] for x in 
                        self.predecessors[node]]).max(axis=0)
                # take inhibitors into account
                if node in self.inhibitors_names:
                    self.values[node] *= 1 - self.inhibitors[node].values

            # take min over all cases 
            count += 1
            
            # check the 2 termination conditions
            #residual = not np.all(abs(outputPrev-newInput))
            self.X1 = pd.DataFrame(self.values)
            #self.X1 = self.values.copy()
            self.debug_values.append(self.X1.copy())
            residual =  (self.X0.fillna(2) - self.X1.fillna(2)).abs().sum().sum()

            #residual = np.nansum(np.square(np.matrix(self.values.values()) - np.matrix(self.X0.values()))) 
            if np.isnan(residual):
                residual = 100
            self.residuals.append(residual)
            #residual = abs(self.X0-self.X1).sum().sum()
        #self.X1 = self.values.copy()
        # set the non-resolved bits to NA
        # TODO
        #newInput[which(abs(outputPrev-newInput) > testVal)] <- NA

    def score(self):
        # with  MMB, fuul bitstrings gives 0.0947746031746032
        # TODO resuse fitness module
        sim = self.X1[self.data.df.columns].values
        data = self.data.df.query("time==@self.time").reset_index(drop=True).values

        #self.sim = sim
        diff = (data-sim)**2
        #self.diff = diff
        S = diff.sum() / float(diff.shape[0] * diff.shape[1])
        return S

    def plotsim(self, time=None):
        import simulator
        sor = simulator.Simulator()

        if time is None:
            time = len(self.debug_values) - 1
        try:
            sor.plot_time_course(self.debug_values[time])
        except:
            sor.plot_time_course(self.X1)

        pylab.tight_layout()



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

