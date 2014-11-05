"""This module is used to simulate multiple runs of Genetic Algorithm using
CellNOptR. It is used to build an average best score vector to be compared with
other type of optimisation (e.g. Metropolis Hasting).


For instance, to run the GA 100 times on the toy model (default value of model
and data parameters)::

    g = BuildAverageGAResults(N=100, model="ToyModelMMB.sif", data="ToyModelMMB.csv")
    g.run()

Default parameters of the GA are used (e.g. popsize=50) but can be changed::

    g.popsize = 20



"""
import numpy
from numpy import linspace, loadtxt, array
#from cellnopt.wrapper import CNORbool
from pylab import *
#from cinapps.mcmc.core import *
#from cinapps.mcmc.diagnostics import *

__all__ = ["BuildAverageGAResults", "MultiGABool", "GABool"]




class GABool(CNORbool):
    def __init__(self, model=None, data=None, verbose=False, popsize=10,
        maxgens=10, mode="T1", debug=True, compression=True, expansion=True):

        #CNORbool.__init__(self, model, data, debug=debug)   
        #self.preprocessing(compression=compression, expansion=expansion)

        self._popsize = popsize
        self._maxgens = maxgens
        self.mode = mode

        self.reset()

    def reset(self):
        self.results = Results(step=self.popsize, N=self.maxgens*self.popsize)
        self.resultsT2 = Results(step=self.popsize, N=self.maxgens*self.popsize)

    def _get_maxgens(self):
        return self._maxgens
    def _set_maxgens(self, N):
        self._maxgens = N
        self.results.N = self._maxgens * self._popsize
    maxgens = property(_get_maxgens, _set_maxgens)

    def _get_popsize(self):
        return self._popsize
    def _set_popsize(self, popsize):
        self._popsize = popsize
        self.results.step = popsize # inside results
        self.results.N = self._maxgens * self._popsize
    popsize = property(_get_popsize, _set_popsize)

    def runT1(self, **kargs):
        self.gaBinaryT1(popsize=self.popsize, maxgens=self.maxgens-1, maxtime=1000000000, **kargs)

        self.results['best_scores'] = self.T1opt.results.Best_score[:]
        self.results['best_score'] = self.T1opt.results.Best_score[-1]
        self.results['best_parameters'] = list(self.T1opt.bString)[:]
        #self.results['all_best_parameters'].append(self.results['best_parameters'])
          

        # same as Best_score but we need to populate scores to prevent plot
        # method to fail.
        self.results['scores'] = list(self.T1opt.results.Best_score_Gen)[:]
        self.results['parameters'] = [[int(y) for y in x.split(',')] for x in  
            self.T1opt.results.Best_bit_Gen]

    def runT2(self, **kargs):
        try:
            self.gaBinaryT2(popsize=self.popsize, maxgens=self.maxgens-1, maxtime=1000000000, **kargs)
        except:
            print "something failed in gaBinaryT2", self.T1opt.results, self.T1opt.bString
        self.resultsT2['best_scores'] = self.T2opt.results.Best_score[:]
        self.resultsT2['best_score'] = self.T2opt.results.Best_score[-1]
        self.resultsT2['best_parameters'] = list(self.T2opt.bString)[:]
        if self.T2opt.results.Best_bit_Gen[0] == "":
            self.resultsT2['parameters'] = [[] for x in self.T2opt.results.Best_bit_Gen]
        else:
            self.resultsT2['parameters'] = [[int(y) for y in x.split(',')] for x in  
                self.T2opt.results.Best_bit_Gen]
                
        #self.resultsT2['all_best_parameters'].append(self.resultsT2['best_parameters'])

        # same as Best_score but we need to populate scores to prevent plot
        # method to fail.
        self.resultsT2['scores'] = list(self.T2opt.results.Best_score_Gen)[:]

    def run(self, mode="T1", **kargs):
        self.runT1(**kargs)
        if mode=="T2": 
            self.runT2(**kargs)

    def plot(self):
        self.results.plot()
  
class MultiGABool(MultiOptimisation):

    def __init__(self, model, data, maxgens=20, popsize=50, Nruns=100, verbose=True):
        super(MultiGABool, self).__init__(model, data, maxgens*popsize, Nruns,
            optim_func=GABool, optim_params={'verbose':verbose,
            'maxgens':maxgens, 'popsize':popsize})

    def run(self, mode="T1", **kargs):
        super(MultiGABool, self).run(mode=mode, **kargs)






class BuildAverageGAResults(object):
    def __init__(self, model="ToyModelMMB.sif", data="ToyModelMMB.csv", N=10,
        popsize=50, maxgens=20, pmutation=0.5, verbose=False):
        """Build an average best score curve over N simulations.

        :param model:
        :param data:
        :param N:
        :param popsize:
        :param maxgens:
        
        To create an instance and run the GA binary 100 times with the default
        GA parameters, type:

            from tools import *
            g = BuildAverageGAResults(N=100, model="ToyModelMMB.sif", data="ToyModelMMB.csv")
            g.run()
            g.plot()

        Results can be saved::

            g.savedata(filename="test.dat")

        And retrieved later on as follows::

            g = BuildAverageGAResults(N=100)
            g.loaddata(filename="test.dat")
            g.plot()
            g.ydata # contains the average best_scores over number of iteration.

        .. note:: The number of iterations is simply popsize times maxgens. There is one
            value for each generation

        Results can be plotted using the R object::

            g.b.plotModel()
            g.b.plotFit()

        """
        assert popsize>2
        assert maxgens>1
        assert pmutation >=0 and pmutation<=1

        self.popsize = popsize
        self.elitism = 5 # default gabinary not to be changed.
        self.pmutation = pmutation

        if self.popsize<= self.elitism:
            self.elitism = self.popsize/2
        self.maxgens = maxgens
        if self.maxgens > 100:
            self.stallgenmax = maxgens
        else:
            self.stallgenmax = 100
        self.N = N

        from cellnopt.data import cnodata
        from cellnopt.wrapper import CNORbool
        self.b = CNORbool(cnodata(model), cnodata(data), verbose=verbose)
        self.data = data
        self.model = model
        self.verbose = verbose
        print "Initialisation done. call run() method and plot() to see the results"
        self.allresults = []
        self.best_bitstrings = []
        self.best_scores = []
        self._computed = False

        self.xdata = [x*self.popsize for x in range(1,self.maxgens+1)]

    def runT1(self, maxtime=100):
        self.run(maxtime=maxtime)

    def run(self, maxtime=100):
        """Run the Genetic Algorithm several times each vector containing the
        scores at each generation.

        """
        if self._computed == True:
            print "Data already computed, set computed to False to force the run()"
        error = 0
      
        for i in range(0, self.N):
            print "Running simulation %s" % i
            # 19 is required instead of 20 because simu stops after the 19th running
            # an extra one anyway making the final value to be 20.
            self.b.run(stallgenmax=self.stallgenmax,
                popsize=self.popsize,maxtime=maxtime, pmutation=self.pmutation,
                maxgens=self.maxgens-1, elitism=self.elitism, verbose=False, show=False, writeall=False)  
            # keep only the first 20 generation (popsize=50) that makes 1000
            # iteration in total.
            res = self.b.T1opt.results.Best_score
          
            # TODO: adapt code to get bewst bitsrint after each run
            self.current_best_bitstring = self.b.T1opt.results.Best_bitString[-1][:]
            self.current_best_bitstring= [int(x) for x in self.current_best_bitstring.split(",")]
         
            if i == 0:
                average = numpy.array(res)
            else:
                try:
                    average += numpy.array(res)
                except:
                    print "Simulation returns a vector of different length. Maybe maxtime is not large enough"
                    # if arrays have different size, += above will not work.
                    # Just increment error to used later when computing the
                    # average. Could happen if maxtime is too short.
                    error+=1
            self.best_scores.append(min(res))
            self.allresults.append(res)
            self.best_bitstrings.append(self.current_best_bitstring)
        self.ydata = average/float(self.N-error)
        self._computed = True
        if len(self.ydata) != len(self.xdata):
            raise ValueError("ydata length does not match xdata length. You may change xdata accordingly or rerun with a larger maxtime parameter.")


    def runT2(self):
        assert self._computed, "must call run() first"



    def plot(self):
        """Plot the scores"""
        from pylab import plot
        plot(self.xdata, self.ydata)
        xlabel("Number of computeScore calls")
        ylabel("Score")
        ylim([0, ylim()[1]])

    def hist(self, nbins=20):
        """Create histogram of the scores"""
        from pylab import hist
        best_scores = [min(x) for x in self.allresults]
        res = hist(best_scores, nbins)
        print "Found %s unique solution:" % len(set(best_scores))
        print set(best_scores)
        return res

    def savedata(self, filename1="BuildAverageGAResults",  filename2="BuildAverageGABitstrings",filename3="BuildAverageGABestScores"):
        """Save the main results (best score, and the averaged bitstring) into a pickle"""
        import pickle
        # b is the CNORbool object that may change in the future so saving it
        # into a pickle is unstable. Let us remove it before saving the pickle
        # and get it back after the pickle is done
        b = self.b
        del self.b
        pickle.dump(self, open(filename1, "w"))
        print "object saved in " + filename1
        self.b = b

        #store best_bitstrings
        pickle.dump(numpy.array(self.best_bitstrings), open(filename2, "w"))
        print "object saved in " + filename2
        
        #store best_scores
        pickle.dump(numpy.array(self.best_scores), open(filename3, "w"))
        print "object saved in " + filename3


    def loaddata(self, filename):
        """load a data set created with savedata method with pickle."""
        import pickle
        res = pickle.load(open(filename))
        #bitstrings = pickle.load(open(filename2))
        self.allresults = res.allresults[:]
        self.xdata = res.xdata[:]
        self.ydata = res.ydata[:]
        self.N = res.N
        self.popsize = res.popsize
        self.maxgens = res.maxgens
        self.computed = True
        data = res.data
        model = res.model
        verbose = res.verbose
        self.verbose = res.verbose
        from cellnopt.data import cnodata
        self.b = CNORbool(cnodata(model), cnodata(data), verbose=verbose)
        del res

    def loadbitstrings(self,filename):
        """load a bitstrings set created with savedata method with pickle."""
        import pickle
       
        bitstrings = pickle.load(open(filename))
        return bitstrings

    def loadbestscores(self,filename):
        """load a best scores from each generation that were created with savedata method with pickle."""
        import pickle

        best_scores = pickle.load(open(filename))
        return best_scores


"""
def get_gabinary_toy("BuildAverageGAResults_toy_20_times_50.dat"):
    g = BuildAverageGAResults(N=100)
    g.loaddata(filename=filename)
    return g.ydata 

def get_gabinary_extliver("BuildAverageGAResults_extliver_20_times_50.dat"):
    g = BuildAverageGAResults(N=100)
    g.loaddata(filename=filename)
    return g.ydata 
"""
