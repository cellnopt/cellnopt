# -*- python -*-
#
#  This file is part of CNO software
#
#  Copyright (c) 2013-2014 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: http://github.com/cellnopt/cellnopt
#
##############################################################################
import os
import sys

import pandas as pd
import pylab
import numpy as np

from cno.core.params import ParamsGA


from cno.core.report import ReportDT
from cno.core.results import DTResults
from cno.core.params import OptionsBase
from cno.core.params import ParamsDT
from cno.core import CNORBase, CNOBase
from biokit.rtools import bool2R

from cno.core.params import params_to_update


__all__ = ["CNORdt"]


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
    def __init__(self, model, data, tag=None, verbose=True,
                 verboseR=False, config=None, use_cnodata=False):
        """

        :param model:
        :param data:
        :param tag:
        :param verbose:
        :param verboseR:
        :param config:
        :return:
        """
        CNOBase.__init__(self,model, data, tag=tag, verbose=verbose,
                        config=config, use_cnodata=use_cnodata)
        self.logging.info("Initialise R session")
        CNORBase.__init__(self, verboseR)

        self._report = ReportDT()
        self._report.Rdependencies = []  # just to speed up code

        self.results = DTResults()

        self.config.General.pknmodel.value = self.pknmodel.filename
        self.config.General.data.value = self.data.filename

        p = ParamsDT()
        self.config.add_section(p)
        self._called = []

    @params_to_update
    def optimise(self, NAFac=1, pmutation=0.5, selpress=1.2, popsize=50,
                 reltol=0.1, elitism=5, maxtime=60, sizefactor=0.0001,
                 time_index_1=1, maxgens=500, maxstallgens=100, bool_updates=10,
                 upper_bound=10, lower_bound=0.8, ga_verbose=True):
        """

        :param lowerB:
        :param upperB:
        :param boolUpdates:

        """

        if int(bool_updates) != len(self.midas.times):
            msg = "boolupdate must be set to number of time points. "
            msg += "Other cases not implemented so far"
            msg += "number time points s %s" % len(self.midas.times)
            print(self.midas.filename)
            raise ValueError(msg)

        # TODO reuse the previous params
        self.logging.info("Running the optimisation. Can take a very long"
                          "time. To see the progression, set verboseR "
                          "attribute to True")
        # update config GA section with user parameters
        self._update_config('DiscreteTime', self.optimise.actual_kwargs)
        self._update_config('GA', self.optimise.actual_kwargs)

        # TODO: add bitstring as input ?

        script = """
        library(CNORdt)
        pknmodel = readSIF("%(pkn)s")
        cnolist = CNOlist("%(midas)s")
        model = preprocessing(cnolist, pknmodel, compression=%(compression)s,
            expansion=%(expansion)s, maxInputsPerGate=%(maxInputsPerGate)s)

        optbs = NULL

        res = gaBinaryDT(CNOlist=cnolist, model=model,
          initBstring=optbs, boolUpdates=%(bool_updates)s,
            popSize=%(popsize)s, maxGens=%(maxgens)s,
            maxTime=%(maxtime)s, elitism=%(elitism)s, pMutation=%(pmutation)s,
            NAFac=%(NAFac)s,  selPress=%(selpress)s, relTol=%(reltol)s, sizeFac=%(sizefactor)s,
            stallGenMax=%(maxstallgens)s, lowerB=%(lower_bound)s,
                  upperB=%(upper_bound)s, verbose=%(ga_verbose)s)

        sim_results = cutAndPlotResultsDT(model=model,
            CNOlist=cnolist, bString=res$bString, boolUpdates=%(bool_updates)s,
            lowerB=%(lower_bound)s, upper=%(upper_bound)s)

        # TODO put this code inside CellNOptR
        signals = colnames(cnolist@signals$`0`)
        # TODO in Dt, MSE are not provided
        # colnames(sim_results$mse) = signals
        #for (i in seq_along(sim_results$simResults[[1]][1,1,])){
        #    colnames(sim_results$simResults[[1]][,,i]) = signals
        #}

        best_bitstring = res$bString
        best_score = res$bScore
        all_scores = res$stringsTolScores
        all_bitstrings = res$stringsTol
        reactions = model$reacID
        results = as.data.frame(res$results)
        stimuli = as.data.frame(cnolist@stimuli)
        inhibitors = as.data.frame(cnolist@inhibitors)
        species = colnames(cnolist@signals[[1]])
        """
        expansion = True
        compression = True

        self._results = {}

        params = {
            'pkn': self.pknmodel.filename,
            'midas': self.data.filename,
             'compression': bool2R(self._compression),
             'expansion': bool2R(self._expansion),
             'cutnonc': bool2R(self._cutnonc),
             'maxInputsPerGate': self._max_inputs_per_gate,

            }

        gad = self.config.GA.as_dict()
        params.update(gad)

        dt_params = self.config.DiscreteTime.as_dict()
        params.update(dt_params)

        params['ga_verbose'] = bool2R(params['ga_verbose'])

        self.session.run(script % params)

        results = self.session.results
        columns_int = ['Generation', 'Stall_Generation']
        columns_float = ['Best_score', 'Avg_Score_Gen', 'Best_score_Gen', 'Iter_time']
        results[columns_int] = results[columns_int].astype(int)
        results[columns_float] = results[columns_float].astype(float)

        reactions = self.session.reactions
        from cno.core.models import DTModels
        try:
            df = pd.DataFrame(self.session.all_bitstrings,
                              columns=reactions)
            self.df = df
        except:
            try:
                df = pd.DataFrame(self.session.all_bitstrings)
                df = df.transpose()
                df.columns = list(reactions)
                self.df = df
            except:
                self.df = pd.DataFrame()

        #df = pd.DataFrame(self.session.all_bitstrings,
        #                      columns=list(self.session.reactions))
        models = DTModels(df)
        models.scores = self.session.all_scores
        models.cnograph.midas = self.data.copy()

        self.best_bitstring = self.session.best_bitstring
        self.species = self.session.species

        results = {
                'best_score': self.session.best_score,
                'best_bitstring': self.session.best_bitstring,
                'all_scores': self.session.all_scores,
                'all_bitstrings': self.session.all_bitstrings,
                'reactions': self.session.reactions,
                'sim_results': self.session.sim_results,  # contains mse and sim at t0,t1,
                'results': results,
                'models': models,
                'stimuli': self.session.stimuli.copy(),
                'inhibitors': self.session.inhibitors.copy(),
                'species': self.session.species,
                #'tag': tag
        }
        results['pkn'] = self.pknmodel
        results['midas'] = self.data

        self.results.results = results
        self.results.models = models

    def plot_errors(self, close=False, show=False):
        results = self.results.results
        midas = results.midas.copy()

        self.simdata = results['sim_results']['simResults'][0]
        self.xCoords = results['sim_results']['optimResults']['xCoords']
        self.yInter = results['sim_results']['optimResults']['yInter']
        self.mse = results['sim_results']['mse']

        Ntimes = len(self.simdata)  # TODO: should match bool_update
        Nspecies = len(self.midas.df.columns)
        Nexp = len(self.midas.experiments.index)
        print Ntimes, Nspecies, Nexp
        #N = Ntimes * Nexp
        #sim = np.array(self.simdata).transpose().reshape(Ntimes,Nexp, Nspecies)
        sim = self.simdata

        # FIXME Ntimes and Nexp may need to be swapped
        #sim = numpy.array(sim).reshape(Ntimes,Nexp,Nspecies)
        species = list(self.session.species)
        df = pd.concat([pd.DataFrame(x, columns=species) for x in sim])
        df.reset_index(inplace=True, drop=True)
        df['experiment'] = list(self.midas.experiments.index) * Ntimes
        df['time'] = [time for time in self.midas.times for x in range(0, Nexp)]
        df['cellLine'] = [self.midas.cellLines[0]] * Nexp * Ntimes

        df.set_index(['cellLine', 'experiment', 'time'], inplace=True)
        df.sortlevel(1, inplace=True)
        #self.df.columns =

        midas.sim = df.copy()
        midas.sim.columns = midas.df.columns
        midas.cmap_scale = 1   # same a CellNOptR
        # need to cut the times

        #valid_times = midas.sim.index.levels[2].values
        #midas.remove_times([x for x in midas.times if x not in valid_times])
        try:midas.plot(mode="mse")
        except:pass

        return midas

    def info(self):
        str_ = "Best bitstring: %s (rmse=%s) " % (self.best_bitstring,self.best_score)
        print(str_)

    def create_report_images(self):

        self.logging.info("Creating figures")
        self._pknmodel.plot(filename=self._report._make_filename("pknmodel.svg"),
                            show=False)
        self._pknmodel.plot(filename=self._report._make_filename("pknmodel.png"),
                            show=False)

        model = self.cnograph.copy()
        model.preprocessing()
        model.plot(filename=self._report._make_filename("expmodel.png"), show=False)

        self.plot_optimised_model(filename=self._report._make_filename("optimised_model.png"),
                                  show=False)

        self.logging.info("Creating mapback model")
        self.plot_mapback_model()
        self._report.savefig("optimised_model_mapback.png")

        self.logging.info("Creating error figure")
        self.plot_errors(show=False)
        self._report.savefig("Errors.png")

        self.logging.info("Creating midas figure")
        self.midas.plot()
        self._report.savefig("midas.png")
        pylab.close()

        self.plot_fitness(show=False, save=True)
        # self._report.savefig("plot_fit.png")

    def create_report(self):

        self._create_report_header()

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
        self._report.add_section(self._report.get_html_reproduce(), "Reproducibility")

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
            txt += "<tr><td>%s:</td><td>%s</td></tr>\n" % (k,v)
        txt += "</table>\n"
        txt += """<img id="img" onclick='changeImage();' src="fit_over_time.png">\n"""
        self._report.add_section(txt, "stats")
        # dependencies
        self._report.write("index.html")

    def cleanup(self):
        # remove the report
        pass

    def _get_stats(self):
        res = {}
        #res['Computation time'] = self.total_time
        try:
            res['Best Score'] = self.bScore
        except:
            pass
        return res

    def _set_simulation(self):

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

    def _init(self):
        print("ARE WE HERE")
        params =  {'pknmodel': self.pknmodel.filename,
                'midas': self.data.filename}
        script = """
        library(CNORdt)
        pknmodel = readSIF("%(pknmodel)s")
        cnolist = CNOlist("%(midas)s")
        model = preprocessing(cnolist, pknmodel)
        """
        script = script % params
        self.session.run(script)

    def simulate(self, bstring, NAFac=1, sizeFac=0.0001, lower_bound=0.8, upper_bound=10,
                 bool_updates=10):
        # given the best bitstring, simulate the data and plot the fit.
        params = {'NAFac': NAFac, 'bool_updates':bool_updates,
                  'sizeFac': sizeFac, 'upper_bound':upper_bound,
                  'lower_bound':lower_bound}
        verboseR = self.verboseR
        self.verboseR = False
        self.session.bstring = bstring
        script = """
        simList = NULL
        indexList = NULL
        mse = computeScoreDT(cnolist, model, bstring, simList, indexList,
            sizeFac=%(sizeFac)s, NAFac=%(NAFac)s, %(bool_updates)s,
            lowerB=%(lower_bound)s, upperB=%(upper_bound)s)
        """
        self.session.run(script % params)
        mse = self.session.mse
        self.verboseR = verboseR
        return mse


def standalone(args=None):
    """This function is used by the standalone application called cellnopt_boolean

    ::

        cellnopt_boolean --help

    """
    if args is None:
        args = sys.argv[:]

    from cno.core.standalone import Standalone
    user_options = OptionsDT()
    stander = Standalone(args, user_options)

    options = stander.options

    if options.onweb is True or options.report is True:
        trainer = CNORdt(options.pknmodel, options.data, verbose=options.verbose,
            verboseR=options.verboseR, config=options.config_file, use_cnodata=options.cnodata)
        trainer.preprocessing()
    else:
        stander.help()

    params = stander.user_options.config.GA.as_dict()
    params.update(stander.user_options.config.DiscreteTime.as_dict())

    trainer.optimise(**params)

    stander.trainer = trainer
    stander.report()


class OptionsDT(OptionsBase):
    def __init__(self):
        prog = "cno_dt"
        version = prog + " v1.0 (Thomas Cokelaer @2014)"
        super(OptionsDT, self).__init__(version=version, prog=prog)
        self.add_section(ParamsGA())
        self.add_section(ParamsDT())


if __name__ == "__main__":
    standalone()
