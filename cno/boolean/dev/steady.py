from cno.core.base import CNOBase
import pandas as pd
import numpy as np

class Steady(CNOBase):

    def __init__(self, model, data, verbose=True):
        super(Steady, self).__init__(model, data, verbose)

        self.time = self.data.times[1]
        self.init(self.time)

    def init(self, time):
        print(time)
        # this should have the same size as experiments ?
        # FIXME is it correct
        N = self.data.df.query('time==@time').shape[0]
        print(N)


        # just a reference to the conditions
        self.inhibitors = self.data.inhibitors
        self.stimuli = self.data.stimuli

        self.inhibitors_names = self.data.inhibitors.columns
        self.stimuli_names = self.data.stimuli.columns


        self.values = {}
        for node in self.model.nodes():
            self.values[node] = np.array([np.nan for x in range(0,N)])

        self.predecessors = {}
        for node in self.model.nodes():
            self.predecessors[node] = self.model.predecessors(node)

        self.successors = {}
        for node in self.model.nodes():
            self.successors[node] = self.model.successors(node)

        self.toflip = [x for x in self.model.edges(data=True) if x[2]['link'] == '-']

        for this in self.stimuli_names:
            self.values[this] = self.stimuli[this].values
        for this in self.inhibitors_names:
            self.values[this] = 1. - self.inhibitors[this].values

        self.and_gates = [x for x in self.model.nodes() if "^" in x]


        self.results = self.data.df.copy()
        self.results = self.results.query("time==@self.time")

    def preprocessing(self):
        self.model.midas = self.data
        self.model.preprocessing()
        self.init(self.time)

    def simulate(self, tick, debug=False):
        
        # First, set the values of the stimuli and inhibitors
        for this in self.stimuli_names:
            self.values[this] = self.stimuli[this].values
        for this in self.inhibitors_names:
            #self.values[this] = self.inhibitors_names[this].values
            self.values[this] = 1. - self.inhibitors[this].values

        self.tochange = [x for x in self.model.nodes() if x not in self.stimuli_names 
                    and x not in self.inhibitors_names and x not in self.and_gates ]

        # what about a species that is both inhibited and measured
        count = 0
        testVal = 1e-3
        nSp = len(self.data.df.columns)
        residual=0

        # FIXME: use numpy.matrix instead of dataframe

        while count < nSp * 1.2 and residual<testVal: 
            self.X0 = pd.DataFrame(self.values)
            # compute AND gates first
            for node in self.and_gates:
                self.values[node] = np.array([self.values[x] for x in 
                        self.predecessors[node]]).min(axis=0)
            for node in self.tochange:
                if len(self.predecessors[node]) == 1:
                    self.values[node] = self.values[self.predecessors[node][0]]
                elif len(self.predecessors[node]) == 0:
                    pass
                else:
                    # ORs gates
                    self.values[node] = np.array([self.values[x] for x in 
                        self.predecessors[node]]).max(axis=0)

            # take min over all cases 
            count += 1
            
            # check the 2 termination conditions
            residual = !all(abs(outputPrev-newInput)
            self.X1 = pd.DataFrame(self.values)
            #residual = abs(self.X0-self.X1).sum().sum()
        # set the non-resolved bits to NA
        # TODO
        #newInput[which(abs(outputPrev-newInput) > testVal)] <- NA



    def score(self):
        
        sim = self.X0[self.data.df.columns]
        data = self.data.df.query("time==@self.time").reset_index(drop=True)

        self.sim = sim
        self.data2 = data
        diff = (data-sim)**2
        S = diff.sum().sum()

        N = diff.shape[0] * diff.shape[1]
        S/=float(N)
        return S









