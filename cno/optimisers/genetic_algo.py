# -*- coding: utf-8 -*-
"""Utilities GA algorithm

Created by Thomas Cokelaer <cokelaer@ebi.ac.uk>
Copyright (c) 2011. GPL


"""
__author__ = """\n""".join(['Thomas Cokelaer <cokelaer@ebi.ac.uk'])

__all__ = ['GABinary']

import copy
import numpy
from easydev import Logging
from cno.misc.profiler import do_profile
import time
import pylab


class GA(Logging):
    def __init__(self, length, verbose=True):
        super(GA, self).__init__(level=verbose)
        self.blength = length



class GABinary(GA):
    """GA algorithm

    ::

        ga = GA(length=16) # define the length
        ga.objective_function = your_function
        ga.run()


    """
    def __init__(self, length, verbose=True, **kargs):
        """**Constructor**


        :param kargs: GA parameters.

        """
        super(GABinary, self).__init__(length, verbose=verbose)


        self.initbstring = kargs.get('initbstring', None)
        self.popsize =  kargs.get('popsize', 50)
        self.pmutation= kargs.get('pmutation', 0.5)
        self.selpress = kargs.get('selpress', 1.2) # use to increase spped of loss of diversity (and convergence) if linear ranking is used between 1 and 2 and if nonlinear between [1, popsize 2]??


        self.elitism = kargs.get('elitism', 5)
        self.reltol = kargs.get('reltol', 0.1)

        #conditions to stop the simulation
        self.maxtime = kargs.get('maxtime', 60)
        self.maxgens = kargs.get('maxgens', 500)
        self.stallgenmax = kargs.get('stallgenmax',100)
        self.stallgenmax = kargs.get('maxstallgen',100)


        self.guess = [1] * self.blength

        self.init()


    def init(self):

        # p = rbind(initBstring,round(matrix(runif(bLength*(PopSize-1)),


        # Generate a bunch of random string, one each population size
        N = self.blength * self.popsize
        self.Pop = numpy.round(numpy.random.uniform(0, 1, N)).reshape(int(self.popsize),
                                                          int(self.blength))
        self.Pop = self.Pop.astype('int64')
        # copy the initBstring in the first row
        self.Pop[0,:] = self.guess

        # we will start with the first column
        self.bestbit = self.Pop[0,:]
        self.bestobj = numpy.inf
        self.bestobj = numpy.inf
        self.stop = False
        self.obj = [0] * self.popsize # used to store all results for each pop
        self.stallgen = 0
        self.g = 0
        #self.res<-rbind(c(g,bestobj,toString(bestbit),stallGen,Inf,Inf,toString(bestbit),0),c(g,bestobj,toString(bestbit),stallGen,Inf,Inf,toString(bestbit),0))
        self.results = {"Generation":[],
                        "Best_score":[],
                        "Best_bitString":[],
                        "Stall_Generation":[],
                        "Avg_Score_Gen":[],
                        "Std_Score_Gen":[],
                        "Best_Score_Gen":[],
                        "Best_bit_Gen":[],
                        "Iter_time":[],
                        "MinScoreGen":[],
                        "MaxScoreGen":[],
                        }

        self.popTol = numpy.array([numpy.nan] * self.blength).reshape(1,self.blength)
        self.popTolScores = numpy.array([numpy.nan])

    #@do_profile()
    def run(self, show=False):
        t0 = time.time()
        start_time = t0
        self.stop = False

        import matplotlib.animation as animation
        if show is True:
            fig, ax = pylab.subplots()
            line, = ax.plot([],[], 'o-', color='b')
            ax.set_ylim(0,1)
            ax.grid()

        while self.stop is False:
            # for each bistring in self.Pop, compute the scores for each string in Pop

            # 99% of the time is spent here in the LiverDREAM data set
            self._scores = [self.getObj(x) for x in self.Pop]
            #print self.g, self._scores, min(self._scores)

            self._rankP = numpy.argsort(self._scores)[::-1] # sorting in decreasing order
            # reorder the matrix according to the ranks
            self.Pop = self.Pop[self._rankP]
            self._scores = numpy.array(self._scores)[self._rankP]

            # Fitness assignment: ranking, linear
            self._fitness = 2 - self.selpress + (2 * (self.selpress-1) *
                                                (numpy.linspace(0, self.popsize-1, self.popsize))
                                                /(self.popsize-1))

            #selection:stochastic uniform sampling
            wheel1 = numpy.cumsum(self._fitness/numpy.sum(self._fitness))
            breaks = numpy.random.uniform(0, 1, 1)*1./self.popsize
            # range stops at self.popsize -1 so that it agrees with CNOR
            breaks = [breaks] + [breaks + float(x)/self.popsize for x in range(1,self.popsize)]

            self._sel = [1] * self.popsize
            for i, x in enumerate(breaks):
                self._sel[i] = numpy.where(numpy.array(wheel1) > x)[0][0]

            # intermediate population
            self.Pop2 = self.Pop[self._sel]
            psize2 = self.Pop2.shape[0]
            psize3 = self.popsize - self.elitism

            self.sizes = [psize2, psize3]
            # Recombination: uniform: each bit has a .5 proba of being
            # inherited from each parent. mates is a 2-column matrice with nrow
            # set to psize3
            # the -1 is used to have a range between 0 and psize2-1 so that
            # indices are correctly set between 0 and array size -1
            _v1 = numpy.ceil(numpy.random.uniform(0, 1, psize3) * psize2)-1
            _v2 = numpy.ceil(numpy.random.uniform(0, 1, psize3) * psize2)-1
            # be zware of the int conversion
            self.mates = numpy.concatenate((_v1, _v2)).astype('int64').reshape(psize3, 2)

             #This holds the probability, for each bit, to be inherited
             # from parent 1 (if TRUE) or 2 (if FALSE)
            _data = numpy.array([x for x in numpy.random.uniform(0 , 1, psize3 * self.blength)])
            _data = _data.reshape(psize3, self.blength)
            self.InhBit = _data < 0.5

            self.Pop3par1 = self.Pop2[self.mates[:,0],:]
            self.Pop3par2 = self.Pop2[self.mates[:,1],:]
            self.Pop3 = self.Pop3par2

            #print self.Pop3par1.shape, self.Pop2.shape, self.mates.shape, self.InhBit.shape
            #print self.InhBit
            self.Pop3par1[self.InhBit]

            self.Pop3[self.Pop3par1[self.InhBit]]

            # mutation
            self.MutProba = numpy.random.uniform(0, 1, psize3*self.blength).reshape(psize3, self.blength)
            self.MutProba = (self.MutProba < (self.pmutation/float(self.blength)))
            self.Pop3[self.MutProba]=  1 - self.Pop3[self.MutProba]

            #Compute stats

            t1 = time.time()

            thisGenBest = self._scores[-1]
            thisGenBestBit = self.Pop[-1, :]

            #if is.na(thisGenBest):
            #    thisGenBest = min(self._scores, na.rm=TRUE)
            #    thisGenBestBit = self.Pop[which(self.scores == thisGenBest)[1],]

            # check that we are not stack somewhere by increasing stallgen
            if thisGenBest < self.bestobj:
                    self.bestobj = thisGenBest
                    self.bestbit = thisGenBestBit
                    self.stallgen = 0
            else:
                self.stallgen += 1

            self.g += 1
            # time aliases

            #save results
            self.results["Generation"].append(self.g)
            self.results["Best_score"].append(self.bestobj)
            self.results["Best_bitString"].append(self.bestbit)
            self.results["Stall_Generation"].append(self.stallgen)
            self.results["Avg_Score_Gen"].append(numpy.mean(self._scores))
            self.results["Std_Score_Gen"].append(numpy.std(self._scores))
            self.results["MinScoreGen"].append(numpy.min(self._scores))
            self.results["MaxScoreGen"].append(numpy.max(self._scores))
            self.results["Best_Score_Gen"].append(thisGenBest)
            self.results["Best_bit_Gen"].append(thisGenBestBit)
            self.results["Iter_time"].append(t1 - t0)

            # breaks

            _tolscore = self._scores[-1] * self.reltol
            _tolbs = numpy.where(self._scores < self._scores[-1] + _tolscore)[0]
            self.info("Generation %s: best score=%s (stall generation: %s)" % (self.g, self.bestobj, self.results['Stall_Generation'][-1]))

            if t1 - start_time > self.maxtime:
                self.info('stopping because time %s exceeds maxtime=%s' % (t1-start_time, self.maxtime))
                self.stop = True
            if self.stallgen > self.stallgenmax:
                self.stop = True
                self.info('Stopping because stall gen max reached the limit')
            if self.g > self.maxgens:
                self.stop = True
                self.info('Stopping because generation max reached the limit')
            if self.bestobj == 0:
                self.stop = True
                self.info('Stopping because score is minimal (deviationPen=0)')

            # update
            #Check for bitstrings that are within the tolerance of the best bitstring



            #, _tolbs

            if len(_tolbs) > 0:
                #print self.popTol
                #print self.Pop[_tolbs, :]
                #FIXME: axis=1 here below will be an error in future numpy version.
                self.popTol = numpy.concatenate( (self.popTol, self.Pop[_tolbs,:]), axis=0)
                self.popTolScores = numpy.concatenate( (self.popTolScores, 
                    self._scores[_tolbs]), axis=1)


            if self.elitism > 0:
                self.Pop = numpy.concatenate( (self.Pop3, self.Pop[self.popsize-self.elitism:self.popsize,:]), axis=0)
            else:
                self.Pop = self.Pop3.copy()

            if show:
                line.set_ydata(self.results['Best_score'])
                line.set_xdata(range(0,self.g))
                ax.relim()
                ax.set_xlim(0, self.g+1)
                ax.set_ylim(0, self.results['MaxScoreGen'][0])
                ax.set_title("Best score=%s" % self.results['Best_score'][-1])
                ax.autoscale_view()
                fig.canvas.draw()
                fig.canvas.flush_events()

            t0 = t1 # essential to stop the while loop

        """
                #Try one point crossover
                #xover<-ceiling(runif(PSize3)*(bLength-1))
                #indices<-matrix(1:bLength,nrow=PSize3,ncol=bLength,byrow=TRUE)
                #InhBit<-matrix(rep(FALSE,PSize3*bLength),nrow=PSize3,ncol=bLength)
                #for(i in 1:PSize3){
                #       InhBit[i,]<-indices[i,]<xover[i]
                #       }
        """

        self.popTol = self.popTol[1:,:]
        self.popTolScores = self.popTolScores[1:]

        _tolbs = numpy.where(self.popTolScores < self._scores[-1] + _tolscore)[0]
        self.popTol = self.popTol[_tolbs,:]
        self.popTolScores = self.popTolScores[_tolbs]
        #self.popTolT = numpy.concatenate( (self.popTol, self.popTolScores), axis=1)
        #self.popTolT = numpy.unique(self.popTolT)#,MARGIN=1)
        #if(!is.null(dim(PopTolT))):
        #PopTol<-PopTolT[,1:(dim(PopTolT)[2]-1)]
        # PopTolScores<-PopTolT[,dim(PopTolT)[2]]

        #else:
        #     PopTol<-PopTolT[1:(length(PopTolT)-1)]
        #     PopTolScores<-PopTolT[length(PopTolT)]


    def plotFit(self):
        import pylab
        pylab.figure(1)
        pylab.clf()
        gen = self.results['Generation']
        bestScore = self.results['Best_Score_Gen']
        avgScore = self.results['Avg_Score_Gen']

        pylab.subplot(2, 1, 1)
        pylab.plot(gen, avgScore)
        pylab.ylabel('Average score of generation')
        ylims = pylab.ylim()
        pylab.ylim(0, ylims[1])

        pylab.subplot(2, 1, 2)
        pylab.plot(gen, bestScore)
        ylims = pylab.ylim()
        pylab.ylim(0, ylims[1])
        pylab.ylabel('Best Score')
        pylab.xlabel('Generation')


    def plotFit2(self):
        from pylab import  clf, plot, vlines, hold, plot, xlabel, ylabel
        clf();
        bestScore = self.results['Best_Score_Gen']
        asg = numpy.array(self.results['Avg_Score_Gen'])
        ssg = numpy.array(self.results['Std_Score_Gen'])
        m1 = numpy.array(self.results['MinScoreGen'])
        m2 = numpy.array(self.results['MaxScoreGen'])
        plot(asg);
        gen = numpy.array(self.results['Generation'])-1
        vlines(gen, asg -ssg, asg+ssg)

        #vlines(gen, m1, m2, color='red')

        hold(True)
        plot(gen, bestScore)
        xlabel('Generation')
        ylabel('Best Score compared to the population average and std')





class GADiscrete(GABinary):

        def __init__(self):
            pass


        #Pop3[MutProba]<-(sample.int(dim(paramsList$type2Funs)[1],length(Pop3[MutProba]),replace = TRUE)) - 1
