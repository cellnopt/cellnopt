"""Old code that may be reused"""
import itertools
import math
import numpy
import random

from pylab import *  


__all__ = ["Optimisation", "MultiOptimisation", "Results", "MultiResults"]


class MultiResults(object):
    """A class to store results from simulation.


    The results of an optimisation are quite similar: scores over time, best
    scores, parameters and so on. So, it is convenient to store those results in
    a well defined objects (see :class:`Results`). When running many simulation,
    each of the results are stored in an individual instance of :class:`Results`.
    
    This class ease the storage of multiple instance of Results and provides
    facilities to extract pertinent information such as average score over the
    different simulation and dedicated plotting tools.

    All results are simply appended to a list :attr:`results`. To add an
    element, it must be an instance of :class:`Results` and you must use the
    :meth:`add_results` method.

    

    """
    def __init__(self):
        self._results = []
        self.graph_options = {}
        self.graph_options['show_individual'] = False

    def _get_N(self):
        """ Return the number of objective function calls"""
        return self._results[0].N
    N = property(_get_N)

    def _get_step(self):
        return self._results[0].step
    step = property(_get_step)

    def _get_xdata(self):
        return range(self.step, self.N+self.step, self.step)
    xdata = property(_get_xdata)

    def add_results(self, res):
        assert isinstance(res, Results)
        self.results.append(res)

    def reset(self):
        self._results = []

    def _get_average_best_score(self):
        """Return the average best score"""
        return mean([x['best_score'] for x in self._results])
    average_best_score = property(_get_average_best_score, 
        doc="return average best score (scalar)")

    def _get_average_best_scores(self):
        """ Return the average best scores over the N calls of the objective function"""
        return numpy.mean([x['best_scores'] for x in self.results], axis=0)
    average_best_scores = property(_get_average_best_scores,
        doc="return average best scores (vector)")
    mean = property(_get_average_best_scores,
        doc="return average best scores (vector)")

    def _get_std_best_scores(self):
        """ Return the standard deviation of the best scores over the N calls of the objective function"""
        return  numpy.std([x['best_scores'] for x in self.results], axis=0)
    std_best_scores = property(_get_std_best_scores,
        doc="return std of best scores (scalar)")
    std = property(_get_std_best_scores,
        doc="return std of  best scores (vector)")

    def _get_results(self):
        return self._results
    results = property(_get_results)

    def _get_best_parameters(self):
        """Return a list of the best parameter for each result"""
        res = [x.best_parameters for x in self.results]
        return res
    best_parameters = property(_get_best_parameters)

    def _get_best_scores(self):
        """Return a list of lists with the all best scores over time for all the runs"""
        a = numpy.array([x['best_scores'] for x in self.results])
        return a
    best_scores = property(_get_best_scores)

    def _get_best_score(self):
        """Return the best score ever"""
        a = numpy.array([x['best_scores'][-1] for x in self.results])
        return a
    best_score = property(_get_best_score)

    def _get_scores(self):
        """Return a list of list which has all the scores from each run."""
        a = numpy.array([x['scores'] for x in self.results])
        return a
    scores = property(_get_scores)

    def _get_parameters(self):
        """Return a list of list which has all the parameters from each run."""
        a = numpy.array([x['parameters'] for x in self.results])
        return a
    parameters = property(_get_parameters)

    def _get_min_index(self):
        """Return the index that the lowest score was achieved for first time."""
        a = [argmin(x) for x in self.best_scores]
        return a
    min_index = property(_get_min_index)

    def plotMulti(self, nbars=10, **kargs):
        label = kargs.get('label', None)
        show_individual = kargs.get('show_individual', self.graph_options['show_individual'])
        legend_size = kargs.get('fontsize',12)
        if show_individual==True:
            for x in self.results:
                plot(self.xdata, x['best_scores'], '--', color="gray")
        x = self.results[0]
        m = self.average_best_scores[:]
        s = self.std_best_scores[:]
        plot(self.xdata, m, linewidth=2, label=label)
        #N = len(m)
        Nm = len(self.xdata)
        step = x.step
        N = x.N
        if nbars >= Nm:
            nbars = Nm
            seprange = self.xdata[:]
        else:
            
            seprange = [int(x) for x in linspace(1, Nm-1, nbars)]
            m = m[seprange]
            s = s[seprange]
            seprange = [step * y for y in seprange]
        #fmt=None to prevent to plot again the m data 
        kargs['label'] = None  # otherwise error bars are added to labels.
        errorbar(seprange, m, yerr=s, fmt=None,ecolor="red",
            elinewidth=2, marker='o', capsize=5, mfc='red', **kargs)
        title("Averaged Best Scores over time (%s runs) " % N)
        xlabel('Iterations')
        ylabel('Score')
        xlim(0, N+step)
        if label: legend(loc=1,prop={'size':legend_size})

    def histScore(self, nbins=50):
        figure()
        data = [x[-1] for x in self.best_scores]
        hist(data, nbins)
        maxScore = max(data)*1.1
        xlim([0, maxScore])
        title('Score Histogram')
        xlabel('Scores')

    def histMinIndex(self, nbins=50):
        """Histogram of the indices at which the best scores were found first.

        This is to investigate the convergence"""
        figure()
        hist(self.min_index, nbins)
        xlabel('Indices where Best Scores was found')
        
    def plot(self, save=True, tag="optimisation", nbins=20,
burnin=0.1,nswaps=1, fontsize=12, label=None):
        """Plots (1) distribution of the scores over time
                 (2) histogram of the scores over time
                 (3) Convergence 
                 (4) Distribution of each bit over time

           If to be used with Greedy search, need to use N=len(scores) instead of self.N

        """
        assert burnin <1, "burnin is in percent and must be less than 100"
        print "Best score is: ", min(self.best_score)
        import numpy

        subplot(2,2,1)
        # ipython --pylab in order to hold the plot each time, otherwise i need hold=True in each plot
        self.plotMulti(label=label)
    
        N = self.results[0].N
        subplot(2,2,2)
        t0 = int(burnin*N)
        hist(self.best_score, bins=nbins, label=label)
        m, M = ylim()
        ylim(m, M*1.1)
        #,label = 'Swaps:%s'%nswaps, alpha=0.5)
        title('Best scores Histogram (%s runs)' % self.N)
        xlabel('Scores')
        if label: legend(loc=1,prop={'size':fontsize}) 
       
        subplot(2,2,3)
        b = numpy.array(self.scores)
        plot(numpy.mean(b,axis=0), label=label)
        #,label = 'Swaps:%s'%nswaps )
        title('Averaged scores over time (%s runs)' % self.N)
        xlabel('Iterations')
        ylabel('Score')
        if label: legend(loc=1,prop={'size':fontsize}) 

        if save:savefig("%s.png" % tag)  # hold the figure for the other plot   


class Results(dict):
    def __init__(self, step=1, N=1):
        super(Results, self).__init__()
        self.step = step
        self.init()
        self._N = N

    def _getN(self):    
        return self._N
    N = property(_getN)

    def _get_xdata(self):
        return range(self.step, self.N+self.step, self.step)
    xdata = property(_get_xdata)

    def copy(self):
        #super(Results, self).copy()
        #self[''] = 
        import copy
        return copy.deepcopy(self)

    def init(self):
        # best_score is the best score at the end of the simulation
        # best_scores is the best score at each iteration step
        # best_parameters is the bitstring corresponding to best score
        # parameters is the list of each bitstrings tried at each iteration step
        # scores
        #self.best_score = 1e6
        self['best_score'] = 1e6
        self['best_parameters'] = []
        self['all_best_parameters'] = []
        self['scores'] = []
        self['best_scores'] = []
        self['parameters'] = []
        self['min_index'] = 0
        self['best_stringT2']=[]

    def _get_scores(self):
        return self['scores']
    scores = property(_get_scores)

    def _get_parameters(self):
        return self['parameters']
    parameters = property(_get_parameters)

    def _get_best_score(self):
        return self['best_score']
    best_score = property(_get_best_score)


    def _get_best_scores(self):
        return self['best_scores']
    best_scores = property(_get_best_scores)

    def _get_best_parameters(self):
        """Return the best parameters that must be provided within the
        optimisation method e.g. in MHIndependent."""
        return self['best_parameters']
    best_parameters = property(_get_best_parameters)
      
    def _get_best_stringT2(self):
        return self['best_stringT2']
    best_stringT2 = property(_get_best_stringT2)

    def plotEdges(self):
        edges = numpy.array([sum(x) for x in self['parameters']])
        plot(edges)
        xlabel('Iterations')
        ylabel('Number of edges on')
        axis([0, len(edges), 0, max(edges)])


    def plot(self, save=False, tag="optimisation", nbins=20, label=None,
            burnin=0.1,nswaps=1, maxScore=None,  legend_size=10):
        """plots 4 subplots

        1. scores versus iteration 
        2. histogram of scores (getting rid of first 10%)
        3. best scores over iteration
        4. bitstrings convergence

        """
        assert burnin <1, "burnin is in percent and must be less than 100"
        print "Best score is: ", min(self['best_scores'])
        import numpy

        scores = self['scores']
        best_scores = self['best_scores'][:]
        bitstrings = self['parameters'][:]


        if maxScore == None:
            maxScore = max(mean(scores)*2, max(scores)*1.1)
        N = len(scores)
        #assert N == self.N

        subplot(2,2,1)
        #plot(self.xdata, scores, label = 'Swaps:%s'%nswaps )
        plot(self.xdata, scores, label=label)
        title('Score Distribution')
        xlabel("$N^o$ of objective function calls")
        ylabel('Score')
        axis([0, self.N, 0, maxScore])
        if label:legend(loc=1,prop={'size':legend_size})

        subplot(2,2,2)
        t0 = int(burnin*N)
        hist(scores[t0:], nbins , label=label)
        #hist(scores[t0:], nbins, label = 'Swaps:%s'%nswaps )
        xlim([0, maxScore])
        m,M = ylim()
        ylim(m,M*1.1)
        title('Score Histogram')
        xlabel('Scores')
        if label:legend(loc=1,prop={'size':legend_size})
        
        subplot(2,2,3)
        #plot(self.xdata, best_scores,label = 'Swaps:%s'%nswaps )
        plot(self.xdata, best_scores, label=label)
        axis([0, self.N, 0, maxScore])
        title('Score Convergence')
        xlabel("$N^o$ of objective function calls")
        ylabel('Score')
        if label:legend(loc=1,prop={'size':legend_size})
        
        # plot all bitstrings. If too many, it is not feasible to plot all of
        # them quickly, so we limit their numbers to 100.
        if N>500:
            step = N/500
        else:
            step = 1
        #print step

        subplot(2,2,4)
        cla() # we do not want to overlap images (pcolor) so clean that subplot
        a = numpy.array(bitstrings)[::step]
        axis([0, a.shape[0], 0, len(bitstrings[0])])
        pcolor(a.transpose())

        xlabel("$N^o$ of objective function calls")
        ylabel('Bits of Bitstring')
        title("Bitstring Convergence of swap %s(Red:1,Blue:0)"%nswaps,fontsize=11)

        if save:savefig("%s.png" % tag)  # hold the figure for the other plot
   


class Optimisation(object):
    """This class provides facilities to run optimisation in the context of
    CellNOpt. 

    :param model:
    :param data:
    :param N: number of calls to objective function.
    :param compression:
    :param expansion:
    :param cutnonc:

    """

    def __init__(self, model, data, N=None, verbose=True, compression=True,
            expansion=True, cutnonc=True):
        #super(Optimisation, self).__init__()
        # user arguments
        if isinstance(model, str):
            model = readSIF(model)
        self.model = model

        if isinstance(data, str):
            self.cnolist = CNOlist(data, verbose=self.verbose)
        else:
            self.cnolist = data


        self._N = N
        self.nTF = 7
        self.step = 1 # delta x; useful for GA results.
        self.verbose = verbose 

        # dedicated to size penalty in objetive function
        self.sizeFac = 0.0001
        self.NAFac = 1

        self.cnolist = rpack_CNOR.CNOlist(self.cnolist)
        self.processed = preprocessing(self.cnolist, self.model, 
            compression=compression, expansion=expansion, cutNONC=cutnonc, verbose=self.verbose)

        self.simlist = prep4sim(self.processed)
        self.indexlist = indexFinder(self.cnolist, self.processed, verbose=False)
        self._params_fuzzy = None
        self._simlist_fuzzy = None
        #self.indexlist = None
        # length of the bitstrings
        #self.m = len(self.model.reacID)
        self._m = len(self.processed.rx2('reacID'))
        self.cardinality = 2**self.m # total number of combination

        # maximum number of iterations
        if N == None:
            self._N = 2**self.m
        else:
            self._N = N

        self.init_bs = [1] * self.m 
        self.init_bs_fuzzy = None

        self.hashscores = {}  # used to store results to avoid recomputation. should have a max size.

        self.reset()

    def _getm(self):
        return self._m
    m = property(_getm, "length of the bitstrings corresponding to the model")

    def _setN(self, value):
        if value<0 or value > self.cardinality:
            raise ValueError("N must be positive and less or equal to cardinality (%s). You provided %s" % (self.cardinality, value))
        else:
            self._N = value
    def _getN(self):    
        return self._N
    N = property(_getN, _setN)


    def reset(self):
        self.results = Results(N=self.N, step=1)
        self.results['best_parameters'] =  self.init_bs[:]
        self.resultsT2 = Results(N=self.N, step=1)
        self.resultsT2['best_parameters'] =  self.init_bs[:]
        self._simlist_fuzzy = None
        self._params_fuzzy = None

    def normalise(self, detection=0, saturation=18000):
        self.cnolist = normaliseCNOlist(self.cnolist, detection=detection, 
            saturation=saturation)



    def __str__(self):
        _str = "Best score over %s iterations: %s" % (self.results['best_score'], len(self.results['scores']))
        _str += "\n"
        _str += "Corresponding best bit string is: %s" % (self.results['best_parameters'])
        _str += "\n"
        return _str



class MultiOptimisation(object):
    """


    :param optim_func: a class representing an objectiv function, with a method
        called run(). func has a results dictionary to be stored Resulst.

    """
    def __init__(self, model, data, N=1000, Nruns=10, optim_func=None,
        optim_params={},verbose=True):
        self.model = model
        self.data = data
        self.N = N
        self.Nruns = Nruns
        self.verbose = verbose 
        optim_params['verbose'] = verbose
        if optim_func == None:
            raise ValueError("optim_func must be provided. Must be a class of type Optimisation")
        try:
            self.func = optim_func(self.model, self.data, self.N, **optim_params)
        except:
            try:
                self.func = optim_func(self.model, self.data,  **optim_params)
            except Exception:
                raise Exception

        self.multi = MultiResults()
        self.multiT2 = MultiResults()

    def reset(self):
        self.multi.reset()
        self.multiT2.reset()

    def run(self, *args, **kargs):
        import time
        self.reset()

        for i in range(0, self.Nruns):

            t1 = time.time()
            #if self.verbose:
            print "Running simulation %s over %s " % (i, self.Nruns)
            self.func.reset()
            self.func.run(*args, **kargs)

            # populate MultiResults with each results structure.
            self.multi.add_results(self.func.results)
            # special case T2 for GA.
            try:
                self.multiT2.add_results(self.func.resultsT2)
            except:
                pass # no T2 data to save.


            t2 = time.time()
            if self.verbose:
                print t2-t1

    def plot(self):
        self.multi.plot()
        try:
            self.multiT2.plot()
        except:
            pass

    def plotModel(self):
        """plot the model/data from the latests run"""
        from rpy2.robjects import FloatVector
        from cellnopt.wrapper import plotModel

        plotModel(self.func.processed.rx2('model'), self.func.cnolist, 
            bString=FloatVector(self.multi.best_parameters))


    def save_average(self):
        a = numpy.array(self.multi.best_scores)
        average = numpy.array(self.multi.mean)
        #m = numpy.min(a, axis=0)
        #M = numpy.max(a, axis=0)
        filename = "results_average.npy"
        print("Saving average best score over time in %s" % filename)
        numpy.save(open(filename, "w"), average)




