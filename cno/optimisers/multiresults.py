"""Old code that may be reused"""
import math
import numpy
import random
import pylab
import pandas as pd

__all__ = ["MultiResults"]


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


    ::

        from cno import Steady, MultiResults
        popsize = 50
        mr = MultiResults(interval=popsize)

        for i in range(0,N):
            s = Steady(self.pknmodel, self.data)
            s.optimise(maxtime=1e6, maxgens=maxgens, popsize=popsize, maxstallgen=200)
            mr.add_scores(s.results.results.best_score)
        mr.plot()


    """
    def __init__(self, interval=1):
        self.scores = []
        self.interval = interval

    def add_scores(self, scores):
        self.scores.append(scores)

    def plot(self, nbars=10, **kargs):
        df = pd.DataFrame(self.scores).T
        mu = df.mean(axis=1)
        sigma = df.std(axis=1)
        N = len(sigma)
        #fmt=None to prevent to plot again the m data 
        #kargs['label'] = None  # otherwise error bars are added to labels.

        xdata = pylab.linspace(self.interval, N*self.interval, N)
        pylab.plot(xdata, mu)
        pylab.errorbar(xdata, mu, yerr=sigma, fmt=None,ecolor="red",
            elinewidth=2, marker='o', capsize=5, mfc='red', **kargs)
        pylab.title("Averaged Best Scores over time (%s runs) " % N)
        pylab.xlabel('Iterations')
        pylab.ylabel('Score')
        pylab.grid()
        #xlim(0, N+step)
        #if label: 
        #    pylab.legend(loc=1,prop={'size':legend_size})

    def hist(self, nbins=50):
        data = [x[-1] for x in self.scores]
        pylab.hist(data, nbins)
        pylab.title('Score Histogram')
        pylab.xlabel('Scores')





