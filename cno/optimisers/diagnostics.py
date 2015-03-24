from pylab import *
import numpy

__all__ = ['Diagnostics']




def runningMovingAverage(data, n=2):
    """Running moving average on a vector
        >>> data = [1,1,3,1,1]
        >>> gen = runningMovingAverage(data)
        >>> list(gen)
        [1.0, 1.0, 2.0, 2.0, 1.0]
    """
    N = len(data)
    NN = float(N)
    for i in range(0,n): yield numpy.mean(data[0:i+1])
    for i in range(n, N): yield numpy.mean(data[i-n+1:i+1])




class Diagnostics(object):
    """Dedicated tools to run diagnostics on optimisation classes

    """

    def __init__(self, *args, **kargs):
        self.reset()
 
    def reset(self):
        self.acceptance = []
        self.alpha = []
        self.acceptanceT2 = []
        self.alphaT2 = [] 

    def _get_acceptanceRate(self):
        return mean(self.acceptance)
    acceptanceRate = property(_get_acceptanceRate, "return mean of list acceptance (acceptance rate)")

    def plotDiagnostics(self, mode="T1"):
        figure(1)
        clf()
        self.plotAcceptance(mode=mode)
        figure(2)
        clf()
        self.plotAutocor()
        figure(3)
        clf()
        self.plotAlpha()

    def plotAlpha(self, hold=False, bins=20):
        if hold == False: 
            clf()
        subplot(2,1,1)
        plot(self.alpha)
        xlabel("iteration")
        ylabel("alpha")
        subplot(2,1,2)
        hist(self.alpha, bins=bins)
        xlabel("alpha")
        ylabel("\#")
        
        if len(self.alphaT2)!=0:
            figure()
            subplot(2,1,1)
            plot(self.alphaT2)
            xlabel("iteration")
            ylabel("alpha")
            subplot(2,1,2)
            hist(self.alphaT2, bins=bins)
            xlabel("alpha")
            ylabel("\#")

    def plotAcceptance(self, lag=None, mode="T1"):
        if mode=="T1":
            if lag == None:
                lag = len(self.acceptance)
            plot(list(runningMovingAverage(self.acceptance, n=lag)))
            xlabel("Iteration")
            ylabel("Running average")
            title("Running Average of the acceptance rate (lag %s)" % lag)
            axis([0, len(self.acceptance), 0, 1])
        else:
            figure()
            if lag == None:
                lag = len(self.acceptanceT2)
            plot(list(runningMovingAverage(self.acceptanceT2, n=lag)))
            xlabel("Iteration")
            ylabel("Running average")
            title("Running Average of the acceptance rate (lag %s)" % lag)
            axis([0, len(self.acceptanceT2), 0, 1])

    def plotAutocor(self, maxlags=None):
        try:
            N = self.N
        except:
            raise ValueError("atribute N not defined.")

        try:
            scores = self.results['scores'][:]
        except:
            raise ValueError("atribute results['scores'] not defined.")

        if maxlags == None:
            maxlags = 50
        acorr(scores[int(N*0.1):N], maxlags=maxlags)
        xlabel("Maxlags")
        ylabel("Autocorrelation")
        title("Model Autocorrelation")

