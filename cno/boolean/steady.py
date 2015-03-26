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



class Steady(CNOBase):
    """Naive implementation of Steady state to help in 
    designing the API"""
    
    def __init__(self, pknmodel, data, verbose=True):
        super(Steady, self).__init__(pknmodel, data, verbose)

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


        self.results = BooleanResults()  # for time T1


        #self.results = self.data.df.copy()
        #self.results = self.results.query("time==@self.time")

        # ignore data of the edge [0:2]
        self.toflip = [x[0:2] for x in self.model.edges(data=True) if x[2]['link'] == '-']

        self.init(self.time)

        
        self.buffering = True

        self.buffer = {}
        self.simulated = {}
        self.debug = True
        self.counter = 0
        self.length_buffer = 10000

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

        for this in self.inhibitors_names:
            self.values[this] = 1. - self.inhibitors[this].values.copy()

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

    def preprocessing(self, expansion=True, compression=True, cutnonc=True):
        self.model.midas = self.data
        self.model.preprocessing(expansion=expansion, compression=compression,
                cutnonc=cutnonc)
        self.init(self.time)
        self.model.buffer_reactions = self.model.reactions
        # we should be using _model from the beginning?
        self._model = self.model

    #@do_profile()
    def simulate(self, tick=1, reactions=None):
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
                # THere are always predecessors
                if length_predecessors[node] != 0:
                    values[node] = bn.nanmin(np.array([values[x] for x in predecessors[node]]), axis=0)
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
                    dummy = np.array([values[x] if (x,node) not in self.toflip 
                        else 1 - values[x] for x in  predecessors[node]])
                    values[node] = bn.nanmax(dummy,  axis=0)

                # take inhibitors into account
                if node in self.inhibitors_names:
                    # if inhibitors is on (1), multiply by 0
                    # if inhibitors is not active, (0), does nothing.
                    values[node] *= 1 - self.inhibitors[node].values
            # here NAs are set automatically to zero because of the int16 cast
            # but it helps speeding up a bit the code by removig needs to take care
            # of NAs. if we use sumna, na are ignored even when 1 is compared to NA            
            self.m1 = np.array([self.previous[k] for k in keys ], dtype=np.int16)
            self.m2 = np.array([values[k] for k in keys ], dtype=np.int16)
            #residual = bn.nansum(np.square(self.m1 - self.m2))
            #residual = np.nansum(np.square(self.m1 - self.m2))
            residual = np.nansum(np.square(self.m1 - self.m2))


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
        # TODO TODO TODO
        #newInput[which(abs(outputPrev-newInput) > testVal)] <- NA
        # loops are handle diffenty

    #@do_profile()
    def score(self, NAFac=1, sizeFac=1e-4):
        # We need also to include NAFac, number of reactions in the model
        # for the sizeFac

        # time 1 only is taken into account
        #self.diff = np.square(self.measures[self.time] - self.simulated[self.time])
        diff = self.measures[self.time] - self.simulated[self.time]
        diff *= diff

        N = diff.shape[0] * diff.shape[1]
        Nna = np.isnan(diff).sum()
        N-= Nna

        #nInTot = number of edges on in global model
        #nInTot = len(self.model.reactions)
        nInTot = self.nInputs # should be correct
        nDataPts = diff.shape[0] * diff.shape[1]
        nDataP = N # N points excluding the NA if any
        #print(N)

        #NAPen = NAFac * sum(self.simulated.isnull())

        # nInTot: number of inputs of expanded miodel
        # nInputs: number of inputs of cut model

        # In CNO:
        # nDataPts = number of points irrespective of NA
        # nDataP  sum(!is.na(CNOlist@signals[[timeIndex]]))
        # nInputs = number of inputs of the cut model

        # for now, ketassume it is the same as the number of reactions
        # TODO AND gates should count for 1 edge
        nInputs = self.number_edges

        sizePen = nDataPts * sizeFac * nInputs / float(nInTot)
        #self.debug("nDataPts=%s" % nDataPts)
        #self.debug("nInputs=%s" % nInputs)
        #self.debug("nInTot=%s" % nInTot)
        #self.debug('sizePen=%s' %sizePen)

        # TODO
        deviationPen = bn.nansum(diff) / 2. # to be in agreement with CNO but wrong
        self.diff = diff / 2.
        #self.debug("deviationPen=%s"% deviationPen)
        #self.debug("Nna=%s"% Nna)
        #self.debug("nDataP=%s"% nDataP)

        deviationPen  /= float(nDataP)
        #self.debug("deviationPen=%s"% deviationPen)
        S = deviationPen + sizePen / nDataP
        return S

    def get_df(self):
        import pandas as pd
        return pd.DataFrame(self.simulated[self.time], columns=self.data.df.columns)

    def plot(self):
        self.model.plot()

    #@do_profile()
    def test(self, N=100):
        # N = 100, all bits on
        # CellNOptR on LiverDREAM  0.85 seconds. 0.58 in cno
        # CellNOptR on LiverDREAM preprocessed) on:0.75 seconds. 1.42 in cno
        # 0.2574948 in CellNOptR  somtimes, we can reach a score=0.019786867202
        # 0.27

        # CellNOptR on ToyMMB      : 0.13        ; 0.22s in cno
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
        system.time(replicate(100,computeScoreT1(cnolist, pknmodel, rep(58) ) ) )
        """
        t1 = time.time()

        reactions = []
        while len(reactions)==0:
            import random
            threshold = np.random.uniform(0,1,1)
            reactions = [r for r in self.model.reactions if random.uniform(0,1)>threshold]

        self.simulate()
        for i in range(0,N):
            #self.init(self.time)
            self.simulate()
            self.score()
        t2 = time.time()
        print(str(t2-t1) + " seconds")        
        return t2-t1

    def plotsim(self, fontsize=16, experiments=None, vmin=0, vmax=1):
        # This is for all experiments is experiments is None
        cm = pylab.get_cmap('gray')
        pylab.clf()

        if experiments is None: # takes the latest (steady state) of each experiments
            data = pd.DataFrame(self.debug_values[-1]).fillna(0.5)
        else:
            data = [(k, [self.debug_values[i][k][experiments] for i in range(0, len(self.debug_values))]) for k in self.debug_values[0].keys()]
            data = dict(data)
            data = pd.DataFrame(data).fillna(0.5)
            data = data.ix[data.index[::-1]]
        self.dummy = data

        pylab.pcolor(data, cmap=cm, vmin=vmin, vmax=vmax,
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
        #pylab.title("Steady state for all experiments(x-axis)\n\n\n\n")
        pylab.tight_layout()

    def plot_errors(self, columns=None):
        # What do we use here: self.values
        print("Use only time 1..")

        # use eval_func with debug one
        debug = self.debug
        self.debug = True
        buffering = self.buffering
        self.buffering = False
        self.eval_func(self.ga.results['Best_bitString'][-1])
        self.buffering = buffering
        self.debug = debug

        if columns is None:
            columns = self.data.df.columns
        X1 = pd.DataFrame(self.debug_values[-1])[columns].copy()
        X1 = self.get_df()
        N = X1.shape[0]

        X1['time'] = [self.time] * N
        X1['cell'] = [self.data.cellLine] * N
        X1['experiment'] = self.data.experiments.index
        X1.set_index(['cell', 'experiment', 'time'], inplace=True)
        self.data.sim.ix[X1.index] = X1
        self.data.plot(mode='mse')
        print("MSE= %s(caspo/cno with only 1 time)" % self.score())
        print("MSE= %s(cellnoptr with only 1 time)" % str(self.score()/2.))

    def optimise2(self, time=None, verbose=True):
        assert len(self.data.times)>=2, "Must have at least 2 time points in the data"
        time1 = self.time
        if time is None:
            self.time = self.data.times[2]
        else:
            self.time = time
        self.init(time)
        prior = list(self.results.results.best_bitstring) # best is the prior
        self.optimise(prior=prior, verbose=verbose)

        # reset to time1
        self.init(time1)

    def plot_errors2(self):
        pass

    def _optimise_gaevolve(self, freq_stats=1, maxgen=200, cross=0.9, elitism=1, 
            mutation=0.02, popsize=80, prior=[]):
        """Example using the library gaevolve

        This was an attempt at using an external library for the GA but it appears
        thst is is slower than our implementation. Dont know why. Could be
        their implementation. Could be the set of parameters.
        """

        # This function is the evaluation function, we want
        # to give high score to more zero'ed chromosomes
        self.count = 0
        def eval_func_in(x):
            return self.eval_func(x)

        self._ga_stats = {'ave':[], 'min':[], 'max':[]}
        def cbcall(ga):
            #print('raw fitness %s' % ga.bestIndividual().score)
            self._ga_stats['min'].append(ga.getStatistics()['rawMin'])
            self._ga_stats['ave'].append(ga.getStatistics()['rawAve'])
            self._ga_stats['max'].append(ga.getStatistics()['rawMax'])


        from pyevolve import G1DList, GSimpleGA, Consts, Selectors
        genome = G1DList.G1DList(len(self.model.reactions))
        genome.evaluator.set(eval_func_in)
        genome.setParams(rangemin=0, rangemax=1, rounddecimal=1e-8, bestrawscore=0)

        ga = GSimpleGA.GSimpleGA(genome)
        ga.setMinimax(Consts.minimaxType["minimize"])
        ga.setGenerations(maxgen)
        ga.setCrossoverRate(cross)
        ga.setMutationRate(mutation)
        ga.setPopulationSize(popsize)
        ga.setElitismReplacement(elitism)
        ga.stepCallback.set(cbcall)
        ga.selector.set(Selectors.GRouletteWheel)
        ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)
        ga.terminationCriteria.set(GSimpleGA.FitnessStatsCriteria)
        #ga.setSortType(Consts.sortType['raw'])

        ga.evolve(freq_stats=freq_stats)
        ga._stats = self._ga_stats.copy()
        #self.plotHistPopScore()
        #self.plotPopScore(hold=hold)
        #from pyevolve import Interaction
        #Interaction.plotHistPopScore(self.ga.getPopulation())
        return ga

    def exhaustive(self):
        from cno.optimisers.binary_tools import permutations
        # create all 
        scores = []
        sizes = []
        from easydev import progress_bar
        N = len(self.model.reactions)
        pb = progress_bar(2**N)
        for i,this in enumerate(permutations(N)):
            self.simulate(this)
            scores.append(self.score())
            pb.animate(i, 0)
            sizes.append(sum(this))
        #self._fill_results()
        self.scores = scores
        self.sizes = sizes
        return scores

    def parameters2reactions(self, chromosome):
        reactions = [x for c,x in zip(chromosome, self._np_reactions) if c==1]
        return reactions

    def reactions2parameters(self):
        pass

    def prior2parameters(self, reactions_on=[], reactions_off=[]):
        prior = [None] * len(self.model.reactions)
        assert len(set(reactions_on).intersection(set(reactions_off))) == 0,\
            "Error. Found reactions in both lists."
        for this in reactions_on:
            idx = self.model.reactions.index(this)
            prior[idx] = 1
        for this in reactions_off:
            idx = self.model.reactions.index(this)
            prior[idx] = 0
        return prior

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

    def optimise(self, verbose=False, maxgens=500, show=False, reltol=0.1, 
        maxtime=60, prior=[], guess=None, reuse_best=True):
        """Using the CellNOptR-like GA"""
        from cno.optimisers import genetic_algo
        ga = genetic_algo.GABinary(len(self.model.reactions), verbose=verbose, 
                maxgens=maxgens, maxtime=maxtime, reltol=reltol)
        # TODO: check length of init guess
        if reuse_best is True:
            try:
                guess = self.results.results.best_bitstring
            except:
                pass

        if guess is not None:
            ga.guess = guess
            ga.init()

        def eval_func_in(x):
            return self.eval_func(x, prior=prior)

        self.counter = 0
        ga.getObj = eval_func_in
        ga.run(show=show)
        self.ga = ga
        self._fill_results()
        return ga

    def _fill_results(self):
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




