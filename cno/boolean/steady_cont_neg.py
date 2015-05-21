from cno.core.base import CNOBase
from cno.core.results import BooleanResults
from cno.core.models import BooleanModels
from cno.misc.profiler import do_profile
from cno.io.reactions import Reaction


import pandas as pd
import numpy as np
import pylab
import time
import bottleneck as bn
from collections import defaultdict
import collections

from cno.boolean.steady import Steady


class SteadyContNeg(Steady):
    """

    
    """
    
    def __init__(self, pknmodel, data, verbose=True):
        super(SteadyContNeg, self).__init__(pknmodel, data, verbose)
        self.results = BooleanResults()  # for time T1


    def init(self, time):
        super(SteadyContNeg, self).init(time)
        self.input_edges = {}
        for i,r in enumerate(self.model.reactions):
            left, right = r.split("=")
            if right in self.input_edges.keys():
                self.input_edges[right].append(i)
            else:
                self.input_edges[right] = [i]


        N = self.N_reactions
        #value_edges = np.zeros((1,N))
        # create function once for all for each edge
        function_edges = [None] * N
        species_in = [None] * N
        for i,reac in enumerate(self.model.reactions):
            if "^" in reac:
                #from cno import Reaction
                #r = Reaction(reac)
                species = self._reactions[i].get_signed_lhs_species()
                def ands(species):
                    species_plus = species['+'][:]
                    species_minus = species['-'][:]
                    data = []
                    if len(species_plus)>0:
                        #print("+")
                        plus = np.min([self.previous[spec]
                            for spec in species_plus], axis=0)
                        data.append(plus.copy())
                    if len(species_minus)>0:
                        #print("-")
                        minus = np.min([1-self.previous[spec]
                            for spec in species_minus], axis=0)
                        data.append(minus.copy())
                    if len(data) == 1:
                        return data[0]
                    else:
                        return np.min(data, axis=0)
                species_in[i] = species.copy()
                function_edges[i] = lambda i: ands(species_in[i]) * self.parameters_in[i]
            else:
                if reac.startswith("!"):
                    species = reac.split("=")[0][1:]
                    function_edges[i] = lambda i: (-self.values[species_in[i]]) * self.parameters_in[i]
                    species_in[i] = species[:]
                else:
                    species = reac.split("=")[0][:]
                    function_edges[i] = lambda i: self.values[species_in[i]] * self.parameters_in[i]
                    species_in[i] = species[:]

        self.species_in = species_in
        self.function_edges = function_edges

    def _init_values(self, time=False):
        self.values = {}
        for node in self.model.nodes():
            # Do we really want NAs ? probably not. fold changes are by
            # definitiotn 0
            #self.values[node] = np.array([np.nan for x in range(0,self.N)])
            self.values[node] = np.zeros(self.N)

        if time > 0:
            self.stimuli = self._stimuli_all
            self.inhibitors = self._inhibitors_all
            for this in self.stimuli_names:
                self.values[this] = self.stimuli[this].values.copy()

            #for this in self.inhibitors_names:
            #    self.values[this] = 1. - self.inhibitors[this].values.copy()

            #toflip = [x[1] for x in self.toflip]
            #for this in toflip:
            #    self.values[this] = np.ones(10)

        else:
            self.stimuli = self._stimuli_none
            self.inhibitors = self._inhibitors_none




    #@do_profile()
    def simulate(self, reactions=None, time=None, ntic=None):
        if time != None:
            assert time in self.data.times
            self.time = time

        # sometimes the shift required to agree with CNOR are not the same at time0 or time1...
        # time T0
        if len(self.toflip):
            self._init_values(0)
            tmp = self._shift
            self._shift = self._shift0
            self._simulate(reactions=reactions, time=0, ntic=ntic)
            self._shift = tmp
            self.dd0 = self.dd.copy()
        else:
            self.simulated[0] = np.zeros(self.measures[0].shape)

        # time T1
        self._init_values(self.time)
        self._simulate(reactions=reactions, time=None, ntic=ntic)

    #@do_profile()
    def _simulate(self, reactions=None, time=None, ntic=None):
        """

        reactions is a list of parameters between 0 and 1
        """
        if reactions is None:
            reactions = [1] * self.N_reactions
        self.parameters_in = reactions[:]
        self.number_edges = self.N_reactions
        if time is None:
            time = self.time
        # what about a species that is both inhibited and measured
        testVal = 1e-3
        import copy
        #values = copy.deepcopy(self.values)
        values = self.values  # !! reference but should be reset when calling _init_values / simulate()


        if self.debug:
            self.debug_values = [values.copy()]
        self.residuals = []
        self.penalties = []

        self.count = 0
        self.nSp = len(values)
        residual = 1.

        frac = 1.2
        # all reactions are always on but for now keep this:
        reactions = self.model.buffer_reactions

        # speed up
        keys = sorted(self.values.keys())

        if ntic is None:
            ntic = self.nSp * frac + self._shift
        else: # we want to use the ntic as unique stopping criteria
            testVal = -1


        # CCK81: inhibitors_failed = ['mTOR', 'PI3K'] but has effects on AKT
        # CCK81: repressors['mTOR'] = ['AKT']
        while ((self.count < ntic) and residual > testVal):
            self.previous = values.copy()

            # 1. for each AND edge and for each normal edge, 
            # get inputs, multiply by strength (parameter)
            #edge_values = [self.function_edges[i](i) for i in range(0,N)]
                
            # 2. for each node, compute max of predecessors
            for node in self.input_edges.keys():

                #data = (self.function_edges[i](i) for i in self.input_edges[node])
                # bootleneck is here. in the for loops and function_edges calls
                M = len(self.input_edges[node])
                if M>1:
                    count = 0
                    data = np.zeros((M, self.N))
                    for i in self.input_edges[node]:
                        data[count] = self.function_edges[i](i)
                        count+=1
                    values[node] = np.sum(data, axis=0)
                else:
                    i = self.input_edges[node][0]
                    values[node] = self.function_edges[i](i)

                if node in self.inhibitors_names and node not in self.inhibitors_failed:
                    # if inhibitors is on (1), multiply by 0
                    # if inhibitors is not active, (0), does nothing.
                    #values[self.drug_targets[node]] *= 1-self.inhibitors[node].values
                    mask = self.inhibitors[node].values == 1
                    values[node][mask] *= 0 # could be optimised as well

                    # could be a value between 0 and 1 or even negative if we have neg values 

                    #if node != 'AKT':
                    #    values[node] = np.clip(values[node], 0, 1)


                # an paradoxical effects induced by drugs ?
                #for inh in self.paradoxical.keys():
                #    if node in self.paradoxical[inh]:
                #        values[node][(self.inhibitors[inh]==1).values] = 1
                for inh in self.repressors.keys():
                    if node in self.repressors[inh]:
                        values[node][(self.inhibitors[inh]==1).values] -= 1



                values[node] = np.clip(values[node], -1, 1)



            # here NAs are set automatically to zero because of the int16 cast
            # but it helps speeding up a bit the code by removig needs to take care
            # of NAs. if we use sumna, na are ignored even when 1 is compared to NA            
            self.m1 = np.array([self.previous[k] for k in keys ], dtype=np.int16)
            self.m2 = np.array([values[k] for k in keys ], dtype=np.int16)
            residual = bn.nansum(np.square(self.m1 - self.m2))
            #residual = np.nansum(np.square(self.m1 - self.m2))


            # TODO stop criteria should account for the length of the species to the
            # the node itself so count < nSp should be taken into account whatever is residual.
            #
            if self.debug:
                self.debug_values.append(self.values.copy())
           
            self.residuals.append(residual)
            if self.stopcount :
                if self.count <10:
                    residual+=1
            self.count += 1

        #if self.debug is True:
        #    # add the latest values simulated in the while loop
        #    self.debug_values.append(values.copy())

        #self._values2 = values

        # Need to set undefined values to NAs

        mask = self.m1 != self.m2
        data = np.array([values[k] for k in keys], dtype=float) 
        data[mask] = np.nan
        self.dd = data

        indices = [keys.index(x) for x in self.data.df.columns]

        if time == 0:
            self.simulated[0] = data[indices,:].transpose()
        else:
            self.simulated[self.time] = data[indices,:].transpose()

    #@do_profile()
    def eval_func(self, chromosome, prior=[]):
        """

        :param prior: a list of same length as chromosome made of 0/1/None

        """
        # TODO limnit the buffering ?
        for i, this in enumerate(prior):
            if this is not None:
                chromosome[i] = this
        
        # using string or tuple takes about the same time but faster than a list
        str_chrome = tuple(chromosome)

        if self.buffering and len(self.buffer)<self.length_buffer and str_chrome in self.buffer.keys():
            return self.buffer[str_chrome]
        else:
            # 110 times faster using numpy array instead of a list...
            reactions = [x for c,x in zip(chromosome, self._np_reactions) if c==1]
            reactions = str_chrome
            self.simulate(reactions=reactions)
            score = self.score()
            if self.buffering is True and len(self.buffer)<self.length_buffer:
                self.buffer[str_chrome] = score
            self.counter +=1
            return score

    #@do_profile()
    def optimise(self, verbose=False, guess=None, 
             elitism=5, method='GA', cooling_rate=0.1,
            nswap=3, popsize=50, mut=0.02, maxgens=50, dimension_bits=1):
        """
        s.optimise(method='inspyred', mut=0.02, pop_size=50, maxgens=40)
        """
        
        # TODO: check length of init guess
        from scipy.optimize import minimize
        from scipy.optimize import basinhopping

        def eval_func_in(x):
            return self.eval_func(x, prior=[])

        N = self.N_reactions

        import inspyred
        from random import Random
        rand = Random()
        if method == 'GA':
            ea = inspyred.ec.GA(rand)
        else:
            ea = inspyred.ec.SA(rand)
        self.scores = []
        self.lower = 0
        self.upper = 1
        self.best_score = 1

        def evaluator(candidates, args):

            dimensions = args['dimensions']
            dimension_bits = args['dimension_bits']
            fitness = []
            for cs in enumerate(candidates):
                #if len(params) != dimensions:
                #print(cs, dimensions, dimension_bits)
                if method == 'GA':
                    params = self._binary_to_real(cs[1], dimensions, dimension_bits)
                else:
                    params = cs[1]
                score = eval_func_in(params)
                #score = eval_func_in(cs)
                fitness.append(score)
                # TODO selection pressure could be added here 
            self.scores.append(min(np.array(fitness)))
            if self.scores[-1] < self.best_score:
                self.best_score =  self.scores[-1]
                i = np.argmin(np.array(fitness))
                self.best_individual = candidates[i]
            start = False
            print(self.scores[-1], self.best_score)

            return fitness
        dimensions = N
        #dimension_bits = 4 # preceision is 1/(2**n-1)
        def generator(random, args):
            # For the first generation only.
            #print('HERE')
            dimension_bits = args['dimension_bits']
            #guess = args['guess'][:]
            pop = [random.choice([0, 1]) 
                    for _ in xrange(self.N_reactions * dimension_bits)]
            #if guess is not None:
            #    pop[0] = guess
            self.pop = pop
            if hasattr(self, 'best_individual'):
                pop = self.best_individual
            return pop

        import inspyred.ec
        import inspyred.ec.ec
        ea.terminator = inspyred.ec.terminators.evaluation_termination
        bounder = inspyred.ec.ec.Bounder(self.lower, self.upper)

        self.storage = []

        def my_observer(population, num_generations, num_evaluations, args):
            #best = max(population)
            #print('{0:6} -- {1} : '.format(num_generations, best.fitness))

            #print(population)
            for pop in population:
                self.storage.append(pop.candidate[:])
        ea.observer = my_observer

        ea.selector = inspyred.ec.selectors.uniform_selection
        #ea.selector = inspyred.ec.selectors.rank_selection
        #ea.selector = inspyred.ec.selectors.fitness_proportionate_selection
        #ea.selector = inspyred.ec.selectors.tournament_selection

        ea.replacer = inspyred.ec.replacers.truncation_replacement
        #ea.replacer = inspyred.ec.replacers.plus_replacement
        self.ea = ea

        final = ea.evolve(generator=generator, evaluator=evaluator, 
                maximize=False, bounder=bounder, temperature=100, cooling_rate=cooling_rate,
                max_evaluations=popsize*maxgens,dimensions=len(self.model.reactions),
                dimension_bits = dimension_bits, guess=guess, gaussian_mean=0.5, gaussian_std=0.5,
                pop_size=popsize, mutation_rate=mut, num_elites=elitism)

        self.ea = ea
        self.final = final
        res = None
        self.best_bitstring = max(self.ea.population).candidate

        try:
            self.best_bitstring = self._binary_to_real(self.best_bitstring,
                len(self.model.reactions), dimension_bits)
        except:
            pass
        return res

    def _binary_to_real(self, binary, dimensions, dimension_bits):
        real = []
        #0 and 1 are the lower/upper bounds
        lo = self.lower
        hi = self.upper
        for d in range(0,dimensions):
            b = binary[d*dimension_bits:(d+1)*dimension_bits]
            real_val = float(int(''.join([str(i) for i in b]), 2))
            value = real_val / (2**(dimension_bits)-1) * (hi - lo) + lo
            real.append(value)

        assert len(real) == dimensions
        return real

    def plot_optimised_model(self):

        ss = self.model.copy()
        for i,r in enumerate(ss.reactions):
            edges = ss.reac2edges(r)
            for edge in edges: 
                ss.edge[edge[0]][edge[1]]['label'] = float(int(self.best_bitstring[i]*100))
                ss.edge[edge[0]][edge[1]]['wc'] = max(5*0.1, self.best_bitstring[i])
                ss.edge[edge[0]][edge[1]]['penwidth'] = max(5*0.1,  5*(self.best_bitstring[i]))
        ss.plot(edge_attribute='wc', cmap='copper')
        return ss

