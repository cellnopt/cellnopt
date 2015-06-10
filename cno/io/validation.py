from cno import CNOGraph
from cno import steady

#Steady



class Validation(object):


    def __init__(self, pknmodel, midas, preprocessing=True):

        print("Computing expected best score")

        self.preprocessing = preprocessing

        #self.c = CNOGraph(pknmodel, midas)
        #self.c.swap_edges(20)

        # compute the real one once for all
        self.real =  steady.Steady(pknmodel, midas)
        if preprocessing is True:
            self.real.preprocessing()
        self.real.optimise(verbose=False)
        self.best_score = self.real.results.results.best_score[-1]
 
        self.best_scores = []
        
    def run(self, N=10, nswap=20, verbose=True, maxstallgen=50, maxtime=60):


        self.sim = steady.Steady(self.real.cnograph, self.real.midas)
        # creates the model, preprocessed
        self.sim.preprocessing()

        from easydev import progress_bar
        pb = progress_bar(N)
        for i in xrange(0,N):
            self.sim = steady.Steady(self.real.cnograph, self.real.midas)
            self.sim.preprocessing()
            self.sim.model.swap_edges(nswap)
            self.sim.preprocessing()
            self.sim.optimise(verbose=verbose, reuse_best=False, 
                    maxstallgen=maxstallgen, maxtime=maxtime)
            score = self.sim.results.results.best_score[-1]
            self.best_scores.append(score)
            pb.animate(i+1)

