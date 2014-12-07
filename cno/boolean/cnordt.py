# -*- python -*-
#
#  This file is part of CORDA software
#
#  Copyright (c) 2014 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: http://www.ebi.ac.uk/~cokelaer/XXX
#
##############################################################################
import os

from easydev import Logging, AttrDict

import pandas as pd
import pylab
import numpy as np

#from cnobool import CNObool, BooleanParameters
from cno.core.params import BooleanParameters
from cno import CNOGraph, XMIDAS
from cno.misc.results import DTResults
from cno.core.report import ReportDT
from cno.core import CNORBase, CNOBase
from biokit.rtools import bool2R

__all__ = ["CNORdt"]




class DTParameters(BooleanParameters):
    # THe keys used here have the same caps as in the R code.
    defaults = {
        'boolupdates': 10,
        'lowerB': 0.8,
        'upperB': 10}
            
    def __init__(self):        
        super(DTParameters, self).__init__()
        self.init_gabinary_t1()
        self.init_dt()

    def init_dt(self):

        self.add_parameter("DT", "--bool-updates", "boolupdates", self.defaults['boolupdates'],
                           "TODO")
        self.add_parameter("DT", "--lowerB", "lowerB", self.defaults['lowerB'],
                           "TODO")
        self.add_parameter("DT", "--upperB", "upperB", self.defaults['upperB'],
                           "TODO")
                           
    
class CNORdt(CNORBase, CNOBase):
    """

    Optimise PKN against data using Discrete Time formalism.


    .. todo:: one difficulty is that boolUpdates
        if different from number of time points. 
        will return a simulation, which size is different from
        the midas file: In cnordt, the data is interpolated
        and replaced but in cellnopt.core, this is 
        more tricky. So, we set boolupdates to be of same size as
        the time points. 
        
    .. warning:: boolUpdates must be same as len(times)
  

    """    
    dt_params = DTParameters.defaults
    def __init__(self, model, data, verbose=True, verboseR=False):

        """.. rubric:: constructor

        """
        CNOBase.__init__(self,model, data, verbose=verbose)
        CNORBase.__init__(self, verboseR)
        self.parameters = {} # fill with GA binary parameters
        self._report = ReportDT()
        self._report._init_report()

        self._optimised = False
        
        params = self.dt_params    
        self.config.add_section("DT")
        self.config.add_option("DT", "boolupdates", params["boolupdates"])
        self.config.add_option("DT", "lowerB", params["lowerB"])
        self.config.add_option("DT", "upperB", params["upperB"])
        
    def optimise(self, boolUpdates, lowerB, upperB, **kargs):
        """

        :param lowerB:
        :param upperB:
        :param boolUpdates:

        """
        params = {
                'boolUpdates': boolUpdates,
                'upperB': upperB,
                'lowerB': lowerB
                }
        self.config.DT.boolupdates  = boolUpdates
        self.config.DT.lowerb  = lowerB
        self.config.DT.upperb  = upperB
        
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
        #for x in mapping.keys():
        #    # update config data structure only if user parameter provided
        #    if x in kargs.keys():
        #        self.config.GA[mapping[x]] = kargs[x]
        #    params[x] = self.config.GA[mapping[x]]
        
           
        script = """
        library(CNORdt)
        pknmodel = readSIF("%(pkn)s")
        cnolist = CNOlist("%(midas)s")

        model = preprocessing(cnolist, pknmodel, compression=%(compression)s,
            expansion=%(expansion)s, maxInputsPerGate=3)

        optbs = NULL

        res = gaBinaryDT(CNOlist=cnolist, model=model,
          initBstring=optbs, boolUpdates=%(boolUpdates)s, maxTime=%(maxtime)s, 
                  lowerB=%(lowerB)s,
                  upperB=%(upperB)s, stallGenMax=%(stallgenmax)s, elitism=%(elitism)s, 
                  popSize=%(popsize)s, sizeFac=1e-4)

        best_bitstring = res$bString
        best_score = res$bScore
        all_scores = res$stringsTolScores
        all_bitstrings = res$stringsTol
        results = as.data.frame(res$results)

        reactions = model$reacID
        stimuli = as.data.frame(cnolist@stimuli)
        inhibitors = as.data.frame(cnolist@inhibitors)
        species = colnames(cnolist@signals[[1]])

        """
        params['maxtime'] =  10
        params['elitism'] = 5 
        params['popsize'] = 30 
        params['stallgenmax'] = 100
        params['pkn'] = self.pknmodel.filename
        params['midas'] = self.data.filename
        #params['tag'] = self.data.filenam
        # FIXME should be generic
        compression = True
        expansion = True
        params['compression'] = bool2R(compression)
        params['expansion'] = bool2R(expansion)
        
        self._results = {}
        self.session.run(script % params)

        results = self.session.results
        columns_int = ['Generation', 'Stall_Generation']
        columns_float = ['Best_score', 'Avg_Score_Gen', 'Best_score_Gen', 'Iter_time']
        results[columns_int] = results[columns_int].astype(int)
        results[columns_float] = results[columns_float].astype(float)


        # for book-keeping, we replace params with actual path (we do not need params anymore)
        #params['model'] = 'PKN-pipeline.sif'
        #params['midas'] = 'MD-pipeline.csv'
        #self._script_optim = script % params
        
        # getting results
        #index = 0
        #self.best_score = self.results['gaBinaryT1'][index].Best_score.min()
        #self.total_time = self.results['gaBinaryT1'][index].Iter_time.sum()
        
        #df = self.results['gaBinaryT1'][index].Best_bitString
        #self.best_bitstring = df.ix[df.index[-1]]
        #self.best_bitstring =  [int(x) for x in self.best_bitstring.split(",")]
        
        #self.optimised_bitstring = {}
        #self.optimised_bitstring['T1'] = self.best_bitstring[:]

        from cno.misc.models import DTModels
        df = pd.DataFrame(self.session.all_bitstrings,
                              columns=list(self.session.reactions))
        models = DTModels(df)
        models.scores = self.session.all_scores
        models.cnograph.midas = self.data.copy()

        results = {
                'best_score': self.session.best_score,
                'best_bitstring': self.session.best_bitstring,
                'all_scores': self.session.all_scores,
                'all_bitstrings': self.session.all_bitstrings,
                'reactions': self.session.reactions,
                'reactions': self.session.reactions,
                'results': results,
                'models':models,
                'stimuli':self.session.stimuli.copy(),
                'inhibitors':self.session.inhibitors.copy(),
                'species':self.session.species,
                #'tag': tag
        }
        results['pkn'] = self.pknmodel
        results['midas'] = self.data

        self.results = DTResults()
        self.results.results = results
        self.results.models = models

    def plot_errors(self, close=False, show=False):
        # todo show parameter

        results = self.results.results
        midas = results.midas

        t0 = results['sim_results']['simResults'][0]['t0']
        t1 = results['sim_results']['simResults'][0]['t1']
        mse = results['sim_results']['mse']
        stimuli = results['stimuli']
        inhibitors = results['inhibitors']
        species = results['species']
        tag = results['tag']

        # need to make sure that order in cellnopt is the same as in the midas
        # shoudl be fine
        # FIXME 40 is hardcoded ...
        assert all([all(midas.stimuli.ix[0] == stimuli.ix[0]) for i in range(0,40)])
        assert all([all(midas.inhibitors.ix[0] == inhibitors.ix[0]) for i in range(0,40)])

        t0 = pd.DataFrame(t0, columns=species)
        t0['experiment'] = midas.experiments.index
        t0['time'] = midas.times[0]
        t0['cell'] = midas.cellLines[0]

        t1 = pd.DataFrame(t1, columns=species)
        t1['experiment'] = midas.experiments.index
        t1['time'] = midas.times[1]
        t1['cell'] = midas.cellLines[0]

        df = pd.concat([t0,t1]).set_index(['cell', 'experiment', 'time'])
        df.sortlevel(1, inplace=True)

        midas.sim = df.copy()
        midas.cmap_scale = 1   # same a CellNOptR
        # need to cut the times

        valid_times = midas.sim.index.levels[2].values
        midas.remove_times([x for x in midas.times if x not in valid_times])
        try:midas.plot(mode="mse")
        except:pass

        return midas

















    def info(self):
        str_ = "Best bitstring: %s (rmse=%s) " % (self.best_bitstring,self.best_score)
        print(str_)

    def create_report_images(self):

        #if self._optimised == False:
        #    raise ValueError("You must run the optimise method first")

        # ust a simple example of settinh the uniprot url
        # should be part of cellnopt.core
        for node in self._pknmodel.nodes():
            self._pknmodel.node[node]['URL'] = "http://www.uniprot.org/uniprot/?query=Ras&sort=score"

        model = self.cnograph.copy()
        model.plot(filename=self._report._make_filename("pknmodel.svg"), show=False)
        model.preprocessing()
        model.plot(filename=self._report._make_filename("expmodel.png"), show=False)

        #self.plot_optimised_model(filename=self._make_filename("optimised_model.png"),
        #                          show=False)

        #c3 = self.get_mapback_model()
        #c3.plot(filename=self._make_filename("optimised_model_mapback.png"),
        #                          show=False)
        #c3 = self.get_mapback_model2()
        #c3.plot(filename=self._make_filename("optimised_model_mapback.png"),
        #        edge_attribute="mycolor", cmap="gray_r")

        self.plot_errors(show=False)
        self._report.savefig("Errors.png")

        self.midas.plot()
        self._report.savefig("midas.png")
        pylab.close()

        self.plot_fitness(show=False, save=True)
        # self._report.savefig("plot_fit.png")

    def plot_fitness(self, show=False, save=True):
        # TODO show and save parameters
        self.results.plot_fit()

        if save is True:
            self._report.savefig("fitness.png")

        if show is False:
            pylab.close()

    def _create_report(self):
        self._report._init_report()
        self._report.directory = self._report.report_directory
        # Save filenames and report in a section
        fname = self._report.directory + os.sep + "PKN-pipeline.sif"

        self.cnograph.to_sif(fname)
        fname = self._report.directory + os.sep + "MD-pipeline.csv"
        self.midas.to_midas(fname)
        txt = '<ul><li><a href="PKN-pipeline.sif">input model (PKN)</a></li>'
        txt += '<li><a href="MD-pipeline.csv">input data (MIDAS)</a></li>'
        txt += '<li><a href="config.ini">Config file</a></li>'
        txt += '<li><a href="rerun.py">Script</a></li></ul>'
        txt += "<bold>some basic stats about the pkn and data e.g. number of species ? or in the pkn section?</bold>"
        self._report.add_section(txt, "Input data files")
        self._report.add_section(
        """
         <div class="section" id="Script_used">
         <object height=120 width=300 type='text/x-scriptlet' border=1
         data="description.html"></object>
         </div>""", "Description")
        txt = """<pre class="literal-block">\n"""
        #txt += "\n".join([x for x in self._script_optim.split("\n") if "write.csv" not in x])
        txt += "o.report()\n</pre>\n"
        self._report.add_section(txt, "Script used")

        txt = """<a href="http://www.cellnopt.org/">
            <object data="pknmodel.svg" type="image/svg+xml">
            <span>Your browser doesn't support SVG images</span> </object></a>"""
        txt += """<a class="reference external image-reference" href="scripts/exercice_3.py">
<img alt="MIDAS" class="align-right" src="midas.png" /></a>"""

        self._report.add_section(txt, "PKN graph", [("http://www.cellnopt.org", "cnograph")])

        self._report.add_section('<img src="expmodel.png">', "Expanded before optimisation")
        self._report.add_section( """
        <img src="optimised_model.png">
        <img src="optimised_model_mapback.png">

        """, "Optimised model")


        self._report.add_section('<img src="fitness.png">', "Fitness")
        self._report.add_section('<img src="Errors.png">', "Errors")

        self._report.add_section(self.get_html_reproduce(), "Reproducibility")
        fh = open(self._report.report_directory + os.sep + "rerun.py", 'w')
        fh.write("from cellnopt.pipeline import *\n")
        fh.write("CNObool(config=config.ini)\n")
        fh.write("c.gaBinaryT1()\n")
        fh.write("c.report()\n")
        fh.close()

        # some stats
        stats = self._get_stats()
        txt = "<table>\n"
        for k,v in stats.iteritems():
            txt += "<tr><td>%s</td><td>%s</td></tr>\n" % (k,v)
        txt += "</table>\n"
        txt += """<img id="img" onclick='changeImage();' src="fit_over_time.png">\n"""
        self._report.add_section(txt, "stats")
        # dependencies
        self._report.write(self._report.report_directory, "index.html")

    def _get_stats(self):
        res = {}
        #res['Computation time'] = self.total_time
        try:
            res['Best Score'] = self.bScore
        except:
            pass
        return res
       
    def _set_simulation(self):
        self.simulate()
        self.midas.create_random_simulation()
        
        Ntimes = self.config.DT.boolupdates
        Nspecies = len(self.midas.df.columns)
        Nexp = len(self.midas.experiments.index)
        print Ntimes, Nspecies, Nexp
        #N = Ntimes * Nexp
        sim = np.array(self.sim).transpose().reshape(Ntimes,Nexp, Nspecies)

        # FIXME Ntimes and Nexp may need to be swapped
        #sim = numpy.array(sim).reshape(Ntimes,Nexp,Nspecies)
        self.df = pd.concat([pd.DataFrame(x) for x in sim])
        self.df.reset_index(inplace=True, drop=True)
        self.df['experiment'] = list(self.midas.experiments.index) * Ntimes
        self.df['time'] = [time for time in self.midas.times for x in range(0, Nexp)]
        self.df['cellLine'] = [self.midas.cellLines[0]] * Nexp * Ntimes
        
        self.df.set_index(['cellLine', 'experiment', 'time'], inplace=True)
        self.df.sortlevel(1, inplace=True)
        #self.df.columns =

        self.midas.sim = self.df.copy()
        self.midas.sim.columns = self.midas.df.columns
        
    def simulate(self):
        """
        
        .. todo:: lowerB, boolupdates should be user parameters as well
        
        
        """
        # given the best bitstring, simulate the data and plot the fit.

        script = """
        library(CNORdt)
        #pknmodel = readSIF("%(pknmodel)s")
        #cnolist = CNOlist("%(midas)s")
        output = cutAndPlotResultsDT(model, %(bs)s, NULL, cnolist, NULL, plotPDF=F,
            boolUpdates=%(boolUpdates)s, lowerB=%(lowerB)s, upperB=%(upperB)s,
            ,sizeFac = 1e-04, NAFac = 1)

        signals = colnames(cnolist@signals[[1]])
        
        #for (i in seq_along(output$simResults[[1]])){
        #    colnames(output$simResults[[1]][[i]]) = signals
        #}

        
        """
        
        params = {
                'boolUpdates': self.config.DT.boolupdates,
                'lowerB': self.config.DT.lowerb,
                'upperB': self.config.DT.upperb,
                'bs': "c(%s)" % ",".join([str(x) for x in self.best_bitstring])}

        script = script % params
        self.session.run(script)
        


