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


class SteadyCont(Steady):
    """

    
    """
    
    def __init__(self, pknmodel, data, verbose=True):
        super(SteadyCont, self).__init__(pknmodel, data, verbose)
        self.results = BooleanResults()  # for time T1

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

        """
        # pandas is very convenient but slower than numpy
        # The dataFrame instanciation is costly as well.
        # For small models, it has a non-negligeable cost.

        # inhibitors will be changed if not ON
        #self.tochange = [x for x in self.model.nodes() if x not in self.stimuli_names
        #            and x not in self.and_gates]

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
        # _shift is set to +1 FIXME +1 is to have same results as in CellnOptR
        # It means that if due to the cycles, you may not end up with same results.
        # this happends if you have cyvles with inhbititions
        # and an odd number of edges. 
        if reactions is None:
            reactions = self.model.buffer_reactions
        self.number_edges = len([r for r in reactions]) + sum([this.count('^') for this in reactions])

        # 10 % time here
        #predecessors = self.reactions_to_predecessors(reactions)
        predecessors = defaultdict(collections.deque)
        for r in reactions:
            k, v = self._reac2pred[r]
            predecessors[k].extend(v)

        # speed up
        keys = sorted(self.values.keys())
        length_predecessors = dict([(node, len(predecessors[node])) for node in keys])

        #self._length_predecessors = length_predecessors
        # if there is an inhibition/drug, the node is 0

        # FIXME is this required ??
        for inh in self.inhibitors_names:
            if length_predecessors[inh] == 0:
                #values[inh] = np.array([np.nan for x in range(0,self.N)])
                values[inh] = np.zeros(self.N)

        # to get same results as in cnor, it is sometimes required
        # to add one more count.
        # to have same results at time 0 as in LiverDream, +3 is required

        if ntic is None:
            ntic = self.nSp * frac + self._shift
        else: # we want to use the ntic as unique stopping criteria
            testVal = -1

        while ((self.count < ntic) and residual > testVal):
            self.previous = values.copy()
            #self.X0 = pd.DataFrame(self.values)
            #self.X0 = self.values.copy()
            # compute AND gates first. why

            # an paradoxical effects induced by drugs ?
            # should be first before updating other nodes
            #for inh in self.paradoxical.keys():
            #    if node in self.paradoxical[inh]:
            #        values[node][(self.inhibitors[inh]==1).values] = 1
            #    #values[inh][(self.inhibitors[inh]==1).values] = 1


            for node in self.and_gates:
                # replace na by large number so that min is unchanged
                # THere are always predecessors
                if length_predecessors[node] != 0:
                    values[node] = bn.nanmin(np.array([self.previous[x] for x in predecessors[node]]), axis=0)
                else:
                    #assert 1==0, "%s %s" % (node, predecessors[node])
                    values[node] = self.previous[node]

            for node in self.tochange:
                # easy one, just the value of predecessors
                #if len(self.predecessors[node]) == 1:
                #    self.values[node] = self.values[self.predecessors[node][0]].copy()
                if length_predecessors[node] == 0:
                    pass # nothing to change
                else:
                    dummy = np.array([self.previous[x] if (x,node) not in self.toflip
                        else 1 - self.previous[x] for x in  predecessors[node]])
                    try:
                        values[node] = bn.nanmax(dummy,  axis=0)
                    except:
                        # in some simple cases, we must reset the type. why.
                        values[node] = bn.nanmax(dummy.astype('int'), axis=0)

                # take inhibitors into account
                if node in self.inhibitors_names and node not in self.inhibitors_failed:
                    # if inhibitors is on (1), multiply by 0
                    # if inhibitors is not active, (0), does nothing.
                    values[node] *= 1 - self.inhibitors[node].values


                # an paradoxical effects induced by drugs ?
                for inh in self.paradoxical.keys():
                    if node in self.paradoxical[inh]:
                        values[node][(self.inhibitors[inh]==1).values] = 1
                for inh in self.repressors.keys():
                    if node in self.repressors[inh]:
                        values[node][(self.inhibitors[inh]==1).values] = 0




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
            self.simulate(reactions=reactions)
            score = self.score()
            if self.buffering is True and len(self.buffer)<self.length_buffer:
                self.buffer[str_chrome] = score
            self.counter +=1
            return score

    #@do_profile()
    def optimise(self, verbose=False, guess=None, method="nelder-mead", elitism=5, 
            nswap=3, pop_size=50, mut=0.02, maxgens=50):
        """Using the CellNOptR-like GA"""
        # TODO: check length of init guess
        from scipy.optimize import minimize
        from scipy.optimize import basinhopping

        from inspyred.ec.utilities import memoize
        def eval_func_in(x):
            #x = [int(x) if x<=1 else 1 for x in x]
            return self.eval_func(x, prior=[])


        if guess is None:
            guess = [1] * len(self.model.reactions)


        def print_fun(x, f, accepted):
            print("at minimum %.4f accepted %d" % (f, int(accepted)))

        N = len(self.model.reactions)

        if method == 'anneal':
            self.da = DiscreteAnnealer(eval_func_in, guess, nswap=nswap)
            state, e = self.da.anneal()
            #state, e = self.da.auto(minutes=1, steps=2)
            print(state, e)
            self.state = state
            self.e = e
            return self.state
        elif method == 'inspyred':
            import inspyred
            from random import Random
            rand = Random()
            ea = inspyred.ec.GA(rand)
            self.scores = []



            #@memoize
            def evaluator(candidates, args):
                fitness = []
                for cs in candidates:
                    score = eval_func_in(cs)
                    fitness.append(score)
                    # TODO selection pressure could be added here 
                self.scores.append(min(np.array(fitness)))
                print(self.scores[-1])
                return fitness

            def generator(random, args):
                # For the first generation only.
                return [random.choice([0, 1]) for _ in range(len(self.model.reactions) * 2)]

            import inspyred.ec
            import inspyred.ec.ec
            ea.terminator = inspyred.ec.terminators.evaluation_termination
            bounder = inspyred.ec.ec.DiscreteBounder([0,1])

            self.generator = generator
            self.evaluator = evaluator

            ea.selector = inspyred.ec.selectors.uniform_selection
            #ea.selector = inspyred.ec.selectors.rank_selection
            #ea.selector = inspyred.ec.selectors.fitness_proportionate_selection
            #ea.selector = inspyred.ec.selectors.tournament_selection

            ea.replacer = inspyred.ec.replacers.truncation_replacement
            #ea.replacer = inspyred.ec.replacers.plus_replacement

            final = ea.evolve(generator=self.generator, evaluator=self.evaluator, 
                    maximize=False, bounder=bounder,
                    max_evaluations=pop_size*maxgens, 
                    pop_size=pop_size, mutation_rate=mut, num_elites=elitism)

            self.ea = ea
            self.final = final
            res = None

        else:
            res = minimize(eval_func_in, guess,
                method=method, bounds=bounds, options={'disp':True})
        
        #self._fill_results()
        return res

    def _fill_results(self):
        # FIXME this could be simplified a lot
        from easydev import AttrDict
        res = AttrDict(**self.ga.results)

        results = pd.DataFrame(self.ga.results)
        columns_int = ['Generation', 'Stall_Generation']
        columns_float = ['Best_score', 'Avg_Score_Gen', 'Best_Score_Gen', 'Iter_time']
        results[columns_int] = results[columns_int].astype(int)
        results[columns_float] = results[columns_float].astype(float)

        results = {
                'best_score': res.Best_score,
                'best_bitstring': res.Best_bitString[-1],
                'all_scores': self.ga.popTolScores,
                'all_bitstrings': self.ga.popTol,
                'reactions': self.model.reactions,
                #'sim_results': self.session.sim_results,  # contains mse and sim at t0,t1,
                'results': results,
                #'models': models,
                #'stimuli': self.session.stimuli.copy(),
                #'inhibitors': self.session.inhibitors.copy(),
                #'species': self.session.species,
        }

        results['pkn'] = self.pknmodel
        results['midas'] = self.data
        #self.results.models = models
        all_bs = self.ga.popTol
        df = pd.DataFrame(all_bs, columns=self.model.reactions)
        models = BooleanModels(df)
        models.scores = results['all_scores']

        self.results.results = results
        self.results.models = models
        self.results.models.cnograph.midas = self.data # to get the MIDAS annotation
        self.best_bitstring = self.results.results.best_bitstring









from simanneal import Annealer
import sys
from simanneal.anneal import time_string

class DiscreteAnnealer(Annealer):
    def __init__(self, eval_func, state, N=1000, nswap=3):
        self.eval_func = eval_func
        super(DiscreteAnnealer, self).__init__(state)
        from easydev import progress_bar as pb
        self.pb = pb(N)
        self.count = 1
        self.steps = N
        self.best_state = self.state[:]
        self.nswap = nswap
        self.best_score = 1

    def move(self):
        #swap states
        # swap 5 states ?
        import random
        self.state = self.best_state[:]
        
        for i in range(0,self.nswap):
            #print self.state
            flipthis = random.randint(0, len(self.state)-1)
            #print flipthis
            self.state[flipthis] = int(self.state[flipthis] + 1) % 2

    def energy(self):
        #self.pb.animate(self.count,0)
        self.count += 1
        score = self.eval_func(self.state)
        if score < self.best_score:
            self.best_score = score
        return score
    
    def update(self, step, T, E, acceptance, improvement):
        elapsed = time.time() - self.start
        if step == 0:
            print(' Temperature        Energy    Accept   Improve     Elapsed   Remaining    Best')
            sys.stdout.write('\r%12.2f  %12.2f                      %s            %s' % \
            (T, E, time_string(elapsed), self.best_score))
            sys.stdout.flush()
        else:
            remain = (self.steps - step) * (elapsed / step)
            sys.stdout.write('\r%12.2f  %12.2f  %7.2f%%  %7.2f%%  %s  %s %s' % \
            (T, E, 100.0 * acceptance, 100.0 * improvement,\
            time_string(elapsed), time_string(remain), self.best_score)),
            sys.stdout.flush()
