import sys
import os

import pandas as pd
import numpy as np
import pylab


#from cnobase import OptionBase
#from cnobool import BooleanParameters
#from htmltools import HTMLReportFuzzy

from cno.boolean import CNOBool   # fuzzy is similar to bool in many aspects


__all__ = ["CNOfuzzy"]



class FuzzyParameters(BooleanParameters):
    # THe keys used here have the same caps as in the R code.
    def __init__(self):        
        super(FuzzyParameters, self).__init__()
        self.init_gabinary_t1()        
    


class CNOfuzzy(CNObool): 
    """


    Optimise first time point in the cnolist and produces a report.



      """
    
    
    gaBinaryT1_params = FuzzyParameters.gaBinaryT1_params

    def __init__(self, pknmodel=None, data=None, config=None,
                 verbose=False, **kargs):
        """.. rubric:: constructor

        """
        super(CNOfuzzy, self).__init__(pknmodel, data, config, formalism="fuzzy",
            verbose=verbose, **kargs)
        
        self.Rdependencies = ["CNORfuzzy"]
        
        self._optimised = False
        if verbose == False:
            self.buffer_Rstdout = True
        self.reset()
 
        self.thresholds = [0.0001, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005,
            0.006, 0.007, 0.008, 0.009, 0.01, 0.013, 0.015, 0.017, 0.02, 0.025, 0.03, 0.05,
            0.1, 0.2, 0.3, 0.5]
 
           
    def reset(self):
        self.results = {
            'gaBinaryT1': []        
        }
            
    def optimise(self, **kargs):
        """alias to :meth:`gaBinaryT1` method"""
        self.gaBinaryT1(**kargs)      
        self._optimised = True

                        
    def gaBinaryT1(self, N=2, **kargs):

        # Save pknmodel in tmp
        # run pipeline
        # read results
        fh1 = self.get_tempfile()
        self.cnograph.export2sif(fh1.name)
        fh1.close()

        # save midas in tmp
        fh2 = self.get_tempfile()
        self.midas.save2midas(fh2.name)
        fh2.close()

        fh_species = self.get_tempfile(suffix=".csv")
        fh_reac = self.get_tempfile(suffix=".csv")
        fh_results = self.get_tempfile()

        params = {
                'pknmodel': fh1.name,
                'midas': fh2.name,
                'fh_species': fh_species.name,
                'fh_reac': fh_reac.name,
                'fh_results': fh_results.name,
                'N':N
                }

        # update config with uesr parameters if provided; keys are user parameter
        # values are internal names used in the config file
        mapping = {
            "selpress": "selection-pressure",
            'reltol': "relative-tolerance",
            'maxtime': "max-time",
            'sizefac': "size-factor",
            'nafac': "na-factor",
            "elitism": 'elitism',
            'popsize': "population-size",
            'stallgenmax': "max-stall-generation",
            'maxgens': "max-generation",
            'pmutation': "probability-mutation",
            'timeindex': "time-index",
            "verbose": "verbose"
        }
        for x in mapping.keys():
            # update config data structure only if user parameter provided
            if x in kargs.keys():
                self.config.GA[mapping[x]] = kargs[x]
            params[x] = self.config.GA[mapping[x]]

        # create the bistring made of ones
        params['bitstring'] =  ",".join(["1" for x in self.cnograph.reacID])
        if params['verbose']==True:
            params['verbose']= "T"
        else:
            params['verbose']= "F"



        script = """
        library(CNORfuzzy)
        pknmodel = readSIF("%(pknmodel)s")
        cnolist = CNOlist("%(midas)s")
        write.csv(list(reacID=pknmodel$reacID), "%(fh_reac)s")
        write.csv(list(species=pknmodel$namesSpecies), "%(fh_species)s")

        paramsList = defaultParametersFuzzy(cnolist, pknmodel)
        paramsList$popSize = %(popsize)s
        paramsList$maxTime = %(maxtime)s
        paramsList$maxGens = %(maxgens)s
        paramsList$stallGenMax = %(stallgenmax)s

        paramsList$optimisation$maxtime = 60*5

        res = CNORwrapFuzzy(cnolist, pknmodel, paramsList=paramsList)

        # 
        #    elitism=%(elitism)s, pMutation=%(pmutation)s,
        #    NAFac=%(nafac)s,  selPress=%(selpress)s, relTol=%(reltol)s, sizeFac=%(sizefac)s,
        #    stallGenMax=%(stallgenmax)s)


        N = %(N)s
        allRes = list()
        paramsList$verbose=TRUE
        for (i in 1:N){
            Res = CNORwrapFuzzy(cnolist, pknmodel, paramsList=paramsList)
            allRes[[i]] = Res
        }
        summary = compileMultiRes(allRes,show=FALSE)

        #sim = plotMeanFuzzyFit(0.01, summary$allFinalMSEs, allRes)

        for (i in seq_along(allRes)){
            write.csv(allRes[[i]]$t1opt$results, paste("%(fh_results)s", "_",i,".csv", sep=""))
        }
        save(allRes, file="allRes.RData")
        save(summary, file="summary.RData")

        """ 
        
        script = script % params
        self._script_optim = script
        self.runRscript(script)


        # for book-keeping, we replace params with actual path (we do not need params anymore)
        params['model'] = 'PKN-pipeline.sif'
        params['midas'] = 'MD-pipeline.csv'
        self._script_optim = script % params



        #reading the saved data
        
        for i in range(1,N+1):
            df = pd.read_csv(fh_results.name + "_" + str(i) + ".csv")
            self.results['gaBinaryT1'].append(df)
        self.reacID = list(pd.read_csv(fh_reac.name).reacID)
        self.species = list(pd.read_csv(fh_species.name).species)
        self.results['reacID'] = self.reacID[:]
        self.results['species'] = self.species[:]
        
        # getting results from the first run only
        index = 1
        self.best_score = self.results['gaBinaryT1'][index].Best_score.min()
        self.total_time = self.results['gaBinaryT1'][index].Iter_time.sum()
        
        df = self.results['gaBinaryT1'][1].Best_bitString
        self.best_bitstring =  df.ix[df.index[-1]]
        self.best_bitstring =  [int(x) for x in self.best_bitstring.split(",")]

        self._cleanup()
        self._optimised = True
        self.optimised_bitstring = {}
        self.optimised_bitstring['T1'] = self.best_bitstring[:]
       

   

    def create_report_images(self):
        if self._optimised == False:
            raise ValueError("You must run the optimise method first")

        
        self._pknmodel.plotdot(filename=self._make_filename("pknmodel.svg"), show=False)
        self.cnograph.plotdot(filename=self._make_filename("expmodel.png"), show=False)
        self.plot_optimised_model(filename=self._make_filename("optimised_model.png"), 
                                  show=False)
 
        self.plot_errors(show=False)

        self.midas.plot()
        self.savefig("midas.png")
        pylab.close()

        self.plot_fitness(show=False, save=True)

    def _plot_fitness(self):
        res = self.results['gaBinaryT1']
        pd.concat([this.Best_score for this in res], axis=1).plot(legend=False)
        pylab.xlabel("Generation")

    
    def report(self, filename="index.html", browse=True, force=False,
               skip_create_images=False):
       
        self.directory = self._init_report()
        
        if skip_create_images == False:
            self.create_report_images()

        report = HTMLReportFuzzy(self)
        self._report(report)

        if browse:
            from browse import browse as bs
            bs(report.directory + os.sep + report.filename)
            
        self.save_config_file()
        
      


    def _set_simulation(self):
        self.simulate()
        self.midas.create_random_simulation()
        
        t0 = np.array(self.sim[[x for x in self.sim.columns if x.startswith("t0")]])
        t0 = pd.DataFrame(t0, columns=self.midas.df.columns)
        t0['experiment'] = self.midas.experiments.index
        t0['time'] = self.midas.times[0]
        t0['cellLine'] = self.midas.cellLines[0]
          
        t1 = np.array(self.sim[[x for x in self.sim.columns if x.startswith("t1")]])
        t1 = pd.DataFrame(t1, columns=self.midas.df.columns)
        t1['experiment'] = self.midas.experiments.index
        t1['time'] = self.midas.times[1]
        t1['cellLine'] = self.midas.cellLines[0]

        df = pd.concat([t0,t1]).set_index(['cellLine', 'experiment', 'time'])
        df.sortlevel(1, inplace=True)
        
        self.midas.sim = df.copy()

    def simulate(self, threshold=0.01, plotPDF=True):
        """
        using plotMeanFuzzyFit
        """
        # given the best bitstring, simulate the data and plot the fit.
        script = """
        library(CNORfuzzy)
        load("allRes.RData") #summary
        load("summary.RData") # allres
        sim = plotMeanFuzzyFit(%(threshold)s, summary$allFinalMSEs, allRes)
        write.csv(sim$simResults, "sim.csv")
        """ % {'threshold': threshold}
        
        self.runRscript(script)
        
        self.sim = pd.read_csv("sim.csv")
        os.remove("sim.csv")
        self._cleanup()
 
    #def compileMultiRes(self, results=None):
    #    self.summary = wrapper_fuzzy.compileMultiRes(self.allRes)


    def plotMSE(self, **kwargs):
        """plot MSEs using interpolation of the results provided by the Fuzzy Analysis"""
        import numpy
        dimRow = self.summary.allFinalMSEs.dim[0]
        allFinalMSEs = numpy.matrix(self.summary.allFinalMSEs)
        allFinalNumParams = numpy.matrix(self.summary.allFinalNumParams)
        catExplore = numpy.zeros((dimRow, len(self.thresholds)))
        AbsNumParams = numpy.zeros((dimRow, len(self.thresholds)))
        AbsMSEs = numpy.zeros((dimRow, len(self.thresholds)))

        # interpolation
        for i in range(0,len(self.thresholds)):
            for j in range(0, dimRow):
                currIX = numpy.where(allFinalMSEs[j,]-allFinalMSEs[j,1]<=self.thresholds[i])[1]
                catExplore[j,i] = numpy.max(currIX)

        for i in range(0,dimRow):
            for j in range(0, len(self.thresholds)):
                AbsNumParams[i,j] = allFinalNumParams[i, catExplore[i,j]]
                AbsMSEs[i,j] = allFinalMSEs[i, catExplore[i,j]]

        # final mean MSEs and number of parameters
        self.meanMSEs = numpy.mean(AbsMSEs, axis=0)
        self.meanNPs = numpy.mean(AbsNumParams, axis=0)

        try:

            fontsize = 20
            fig1 = pylab.figure()
            ax1 = fig1.add_subplot(111)
            line1 = ax1.semilogx(self.thresholds, self.meanMSEs, 'b-o', **kwargs)
            pylab.ylabel("MSEs", fontsize=fontsize)
            pylab.xticks(fontsize=16)
            pylab.yticks(fontsize=16)

            pylab.axis([self.thresholds[0], self.thresholds[-1],
                min(self.meanMSEs)/1.01,max(self.meanMSEs)*1.01])

            ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)
            line2 = ax2.plot(self.thresholds, self.meanNPs, 'r-o')
            ax2.yaxis.tick_right()
            ax2.yaxis.set_label_position("right")
            pylab.ylabel("Number of Parameters", fontsize=fontsize)

            pylab.legend((line1, line2), ("mean MSEs", "mean Number of Parameters"),
                loc="center left", prop={'size':13})
            pylab.grid()
            pylab.xticks(fontsize=16)
            pylab.yticks(fontsize=16)
        except ImportError, e:
            print(e)
            print("install pylab to use this function")
        #show()





def standalone(args=None):
    """This function is used by the standalone application called cellnopt_boolean

    ::

        cellnopt_boolean --help

    """
    if args == None:
        args = sys.argv[:]

    user_options = OptionFuzzy(prog="cellnopt_fuzzy")

    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    else:
        options = user_options.parse_args(args[1:])

    o = CNOfuzzy(options.model, options.data, verbose=options.verbose)
    o.optimise()

    if options.report:
        o.report()
    else:
        print("No report request (use --report)")


class OptionFuzzy(OptionBase):

    def  __init__(self, version="1.0", prog=None):
        usage = """usage: python %s --data ToyModelMMB.csv --model ToyModelMMB.sif""" % prog
        super(OptionFuzzy, self).__init__(usage=usage, version=version, prog=prog)
        self.add_gaBinaryT1_options()

    def add_gaBinaryT1_options(self):
        """The input options.

        Default is None. Keep it that way because otherwise, the contents of
        the ini file is overwritten in :class:`apps.Apps`.
        """
        params = FuzzyParameters()
        group = self.add_argument_group("Genetic Algorithm", 
                    """This section gathers the parameters of the Genetic Algorithm
                    """)
        keys = params.get_keys_from_section("GA")
        
        for key in  keys:
            param = params.parameters[key]
            
            kargs = param._get_kargs()
            help = str(kargs['help'])
            del kargs["help"]
            print(help)
            group.add_argument(param.name , help=""+help, **kargs)
        

if __name__ == "__main__":
    """Used by setup.py as an entry point to :func:`standalone`

    """
    standalone(sys.argv)

