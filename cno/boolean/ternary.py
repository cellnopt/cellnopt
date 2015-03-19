from cno.core.base import CNOBase

import pandas as pd
import numpy as np
import pylab
from cno.misc.profiler import do_profile
from cno.io.reactions import Reaction
import time
import bottleneck as bn
from collections import defaultdict
import collections
from cno.boolean.steady import Steady



class Ternary(Steady):
    """Naive implementation of Steady state to help in 
    designing the API"""
    
    def __init__(self, pknmodel, data, verbose=True):
        super(Ternary, self).__init__(pknmodel, data, verbose)

        self.model = self.pknmodel.copy()
        # to speed up code later on
        self.model.buffer_reactions = self.model.reactions[:]
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

        # note that the order of the rows is experiments as defined in
        # data.df not data.experiments
        self.measures = {}
        # FIXME No need for time zero but if so, need to re-order the experiments
        #self.measures[0] = self.data.df.query("time==0").reset_index(drop=True).values
        df = self.data.df.query("time==@self.time")
        df = df.ix[self.data.cellLine]
        df = df.ix[self.stimuli.index]
        df = df.reset_index(drop=True).values
        self.measures[self.time] = df
        
        self.buffering = True

        self.buffer = {}
        self.simulated = {}
        self.debug = True

    def init(self, time):
        # Here, we define the values of the stimuli and inhibitors
        # based on the order provided inside self.stimuli and 
        # self.inhibitors

        # Later, one has to be cautious with the measured data, which 
        # order may be different !!

        assert time in self.data.times
        self.values = {}
        for node in self.model.nodes():
            # Do we really want NAs ? probably not. fold changes are by
            # definitiotn 0
            #self.values[node] = np.array([np.nan for x in range(0,self.N)])
            self.values[node] = np.array([0 for x in range(0,self.N)])

        for this in self.stimuli_names:
            self.values[this] = self.stimuli[this].values.copy()

        # TODO what should be the initial values here ? Let us set the values to -1
        # TODO/TODO
        # TODO/TODO
        # TODO/TODO
        for this in self.inhibitors_names:
            self.values[this] =  1 - self.inhibitors[this].values.copy()
            self.values[this] =   self.inhibitors[this].values.copy()

        self.and_gates = [x for x in self.model.nodes() if "^" in x]

        self.predecessors = {}
        for node in self.model.nodes():
            self.predecessors[node] = self.model.predecessors(node)
        self.number_predecessors = {}
        for node in self.predecessors.keys():
            self.number_predecessors = len(self.predecessors[node])

        self.successors = {}
        for node in self.model.nodes():
            self.successors[node] = self.model.successors(node)

        self.nInputs = np.sum([len(self.model.predecessors(x)) for x in self.model.nodes()])
        self.nInputs -= len(self.model._find_and_nodes()) # get rid of the output edge on AND gates

        self.tochange = [x for x in self.model.nodes() if x not in self.stimuli_names
                    and x not in self.and_gates]

        self._reactions = [Reaction(r) for r in self.model.reactions]

        self._np_reactions = np.array(self.model.reactions)

        self._reac2pred = {}
        for r in self.model.reactions:
            reac = Reaction(r)
            if "^" in reac.lhs:
                self._reac2pred[r] = (r, reac.lhs_species)
            else:
                self._reac2pred[r] = (reac.rhs, reac.lhs_species)

    #@do_profile()
    def simulate(self, tick=1, debug=False, reactions=None):
        """

        """
        # pandas is very convenient but slower than numpy
        # The dataFrame instanciation is costly as well.
        # For small models, it has a non-negligeable cost.

        # inhibitors will be changed if not ON
        #self.tochange = [x for x in self.model.nodes() if x not in self.stimuli_names
        #            and x not in self.and_gates]

        # what about a species that is both inhibited and measured
        testVal = 1e-3

        values = self.values.copy()

        if self.debug:
            self.debug_values = []
        self.residuals = []
        self.penalties = []

        self.count = 0
        self.nSp = len(values)
        residual = 1.

        frac = 1.2
        # #FIXME +1 is to have same resrults as in CellnOptR
        # It means that if due to the cycles, you may not end up with same results.
        # this happends if you have cyvles with inhbititions
        # and an odd number of edges. 
        if reactions is None:
            reactions = self.model.buffer_reactions
        self.number_edges = len(reactions)

        # 10 % time here
        #predecessors = self.reactions_to_predecessors(reactions)
        predecessors = defaultdict(collections.deque)
        for r in reactions:
            k,v = self._reac2pred[r]
            predecessors[k].extend(v)

        # speed up
        keys = self.values.keys()
        length_predecessors = dict([(node, len(predecessors[node])) for node in keys])


        #self._length_predecessors = length_predecessors
        # if there is an inhibition/drug, the node is 0
        values = self.values.copy()
        for inh in self.inhibitors_names:
            if length_predecessors[inh] == 0:
                #values[inh] = np.array([np.nan for x in range(0,self.N)])
                #values[inh] = np.array([0 for x in range(0,self.N)])
                values[inh] = np.zeros(self.N)

        while (self.count < self.nSp * frac +1.) and residual > testVal: 
            self.previous = values.copy()
            #self.X0 = pd.DataFrame(self.values)
            #self.X0 = self.values.copy()
            # compute AND gates first. why
            for node in self.and_gates:
                # replace na by large number so that min is unchanged
                if length_predecessors[node] != 0:
                    # not a min anymore but a consensus
                    #values[node] = np.nanmin(np.array([values[x].copy() for x in
                    #    predecessors[node]]), axis=0)
                    #print(np.array([values[x].copy() for x in
                    #    predecessors[node]]))
                    # TODO/TODO
                    # TODO/TODO
                    # TODO/TODO
                    dummy = np.array([values[x].copy() for x in
                        predecessors[node]])

                    values[node] = consensus2(dummy.transpose())
                else:
                    values[node] = self.previous[node]

            for node in self.tochange:
                # easy one, just the value of predecessors
                #if len(self.predecessors[node]) == 1:
                #    self.values[node] = self.values[self.predecessors[node][0]].copy()
                if length_predecessors[node] == 0:
                    pass # nothing to change
                else:

                    # TODO/TODO
                    # TODO/TODO
                    # here the nagative edges are not 1-x anymore but just -x
                    dummy = np.array([values[x] if (x,node) not in self.toflip
                        else - values[x] for x in  predecessors[node]])
                    # TODO/TODO not an AND but a consensus
                    # TODO/TODO
                    # TODO/TODO
                    values[node] = accept_anything(dummy.transpose())

                # take inhibitors into account
                # TODO/TODO do we want to change this behaviour ?
                # TODO/TODO
                # TODO/TODO
                if node in self.inhibitors_names:
                    values[node] *=  1- self.inhibitors[node].values
            # 30 % of the time is here 
            # here NAs are set automatically to zero because of the int16 cast
            # but it helps speeding up a bit the code by removig needs to take care
            # of NAs. if we use sumna, na are ignored even when 1 is compared to NA            
            self.m1 = np.array([self.previous[k] for k in self.previous.keys() ], dtype=np.int16)
            self.m2 = np.array([values[k] for k in self.previous.keys() ], dtype=np.int16)
            residual = bn.nansum(np.square(self.m1 - self.m2))

            # TODO stop criteria should account for the length of the species to the
            # the node itself so count < nSp should be taken into account whatever is residual.
            #
            if self.debug:
                self.debug_values.append(self.previous.copy())
           
            self.residuals.append(residual)
            self.count += 1

        if self.debug is True:
            # add the latest values simulated in the while loop
            self.debug_values.append(values.copy())

        # Need to set undefined values to NAs
        self.simulated[self.time] = np.array([values[k] 
            for k in self.data.df.columns ], dtype=float)#.transpose()

        self.prev = {}
        self.prev[self.time] = np.array([self.previous[k] 
            for k in self.data.df.columns ], dtype=float)#.transpose()

        mask = self.prev[self.time] != self.simulated[self.time]
        self.simulated[self.time][mask] = np.nan

        self.simulated[self.time] = self.simulated[self.time].transpose()

        # set the non-resolved bits to NA
        # TODO
        #newInput[which(abs(outputPrev-newInput) > testVal)] <- NA
        # loops are handle diffenty

   
   
    def plotsim(self, fontsize=16, experiments=None):
        super(Ternary, self).plotsim(vmin=-1)

   
def accept_anything(x, axis=0):
    res = []
    for this in x:

        m1 = min(this)
        m2 = max(this)
        if m1 >= 0 and m2>=0:
            res.append(max(m1,m2))
        elif m1 <=0 and m2<=0:
            res.append(min(m1,m2))
        else:
            res.append(0)
    res = np.array(res)
    res.reshape(x.shape[0],1)
    return res


def consensus2(x):
    res = []
    for this in x:
        if all(this==1):
            res.append(1)
        elif all(this==-1):
            res.append(-1)
        else:
            res.append(0)
    res = np.array(res)
    res.reshape(x.shape[0],1)
    return res

def consensus_ternary(x):
    # 2.24us
    if len(x)>2:
        return consensus_ternary([consensus_ternary(x[0:2]), x[2:]])
    else:
        if x[0] == x[1]:
            return x[0]
        else:
            return 0
        return consensus

def consensus3(x):
    #using np array as inut 26us
    #using list 401 ns
    # 26us
    if len(x)>2:
        return consensus3([consensus3(x[0:2]), x[2:]])
    else:
        return x[0] & x[1] | ( (x[0]!=-1)&0) | ((x[1]!=-1) &0)
