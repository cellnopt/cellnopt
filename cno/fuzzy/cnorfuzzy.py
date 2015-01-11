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
import sys
import os

import pandas as pd
import numpy as np
import pylab

from cno.core import CNOBase
from cno.core.params import OptionsBase

from cno.core import CNORBase
from cno.core.results import FuzzyResults
from cno.core.report import ReportFuzzy
from cno.core.models import FuzzyModels
from cno.core.params import ParamsGA, ParamsFuzzy

from biokit.rtools import bool2R

from cno.core.params import params_to_update


__all__ = ["CNORfuzzy"]

"""
$redThresh  0e+00 1e-04 5e-04 1e-03 3e-03 5e-03 1e-02

"""

class FuzzyParameters(ParamsGA):
    # THe keys used here have the same caps as in the R code.
    def __init__(self):
        super(FuzzyParameters, self).__init__()
        self.nTF = 7


class CNORfuzzy(CNOBase, CNORBase):
    """

    Optimise first time point in the cnolist and produces a report.

    """
    def __init__(self, model=None, data=None, tag=None, verbose=True,
                 verboseR=False, config=None, use_cnodata=False):
        """.. rubric:: constructor

        :param model:
        :param data:
        :param verbose:
        :param verboseR:
        :param config:
        :return:

        """
        CNOBase.__init__(self,model, data, verbose=verbose, 
                config=config, use_cnodata=use_cnodata)
        self.logging.info("Initialise R session")
        CNORBase.__init__(self, verboseR=verboseR)

        self._report = ReportFuzzy()
        self._report.Rdependencies = []  # just to speed up report

        self.results = FuzzyResults()

        self.config.General.pknmodel.value = self.pknmodel.filename
        self.config.General.data.value = self.data.filename

        p = ParamsFuzzy()
        self.config.add_section(p)

        self.thresholds = [0.0001, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005,
            0.006, 0.007, 0.008, 0.009, 0.01, 0.013, 0.015, 0.017, 0.02, 0.025, 0.03, 0.05,
            0.1, 0.2, 0.3, 0.5]

    def _load_cnolist(self):
        script_template = """
        library(CNORfuzzy)
        pknmodel = readSIF("%(pkn)s")
        cnolist = CNOlist("%(midas)s")"""
        script = script_template % {
                'pkn': self.pknmodel.filename,
                'midas': self.data.filename}
        self.session.run(script)

    @params_to_update
    def optimise(self, N=2,
        NAFac=1, pmutation=0.5, selpress=1.2, popsize=50,
        reltol=0.1, elitism=5, maxtime=60, sizefactor=0.0001,
        time_index_1=1, maxgens=500, maxstallgens=100, ga_verbose=True, **kargs):

        self.logging.info("Running the optimisation. Can take a very long"
                          "time. To see the progression, set verboseR "
                          "attribute to True")
        # update config GA section with user parameters
        self._update_config('GA', self.optimise.actual_kwargs)
        self._update_config('Fuzzy', self.optimise.actual_kwargs)

         # keep track of the GA parameters, which may have been update above
        gad = dict([(k, self.config.GA[k].value)
            for k in self.config.GA._get_names()])

        fuzzyd = dict([(k, self.config.Fuzzy[k].value)
            for k in self.config.Fuzzy._get_names()])

        print(N)
        print(fuzzyd)
        script = """
        library(CNORfuzzy)
        cnolist = CNOlist("%(midas)s")
        pknmodel = readSIF("%(pkn)s")

        # pknmodel is processed internally. Need to change the R API
        #model = preprocessing(cnolist, pknmodel, compression=%(compression)s,
        #    expansion=%(expansion)s, maxInputsPerGate=3)

        # which one to use ?? pknmodel or model
        # looks like CNORwrap overwrite paramsList$model
        # with the processed one
        # $model not used in the function so one can provide anything
        paramsList = defaultParametersFuzzy(cnolist, pknmodel)
        paramsList$popSize = %(popsize)s
        paramsList$maxTime = %(maxtime)s
        paramsList$maxGens = %(maxgens)s
        paramsList$elitism = %(elitism)s
        paramsList$stallGenMax = %(maxstallgens)s
        paramsList$optimisation$maxtime = %(optimisation_max_time)s

        N = %(N)s
        allRes = list()
        paramsList$verbose=TRUE
        for (i in 1:N){
            Res = CNORwrapFuzzy(cnolist, pknmodel, paramsList=paramsList)
            allRes[[i]] = Res
        }
        summary = compileMultiRes(allRes,show=FALSE)
        summary = compileMultiRes(allRes,show=T)
        # signals order is not sorted in CellNOptR
        signals = colnames(cnolist@signals[[1]])
        #sim = plotMeanFuzzyFit(0.01, summary$allFinalMSEs, allRes)

        res1 = allRes[[1]]
        best_score = res1['currBestDiscrete']
        best_bitstring = res1['intString']
        """

        expansion = True
        compression = True

        params = {
            'pkn': self.pknmodel.filename,
            'midas': self.data.filename,
            'compression': bool2R(compression),
            'expansion': bool2R(expansion)
            }
        params.update(gad)
        params.update(fuzzyd)

        print(params)
        self.session.run(script % params)
        allRes = self.session.allRes

        # The contents of allRes is a list of N Res structures
        # Each Res structure is itself made of 9 structures with those keys:
        # 't1opt': ['stringsTol', 'bScore', 'results', 'currBest', 'bString', 'stringsTolScores']
        # 'paramsList': parameters used in particular GA params, optim params, inputs
        # 'unRef':
        # 'processedModel':
        # 'currBestDiscrete': best score
        # 'bit'
        # 'intString'
        # 'redRef'
        # 'cutBit'

        res1 = allRes[0] # same for all indices
        #reactions = res1['paramsList']['model']['reacID']
        species = res1['paramsList']['model']['namesSpecies']

        reactions = res1['processedModel']['reacID']

        # TODO: find best MSE/bitstring amongst the N runs
        # redRef contains a MSE for each threshold
        res1['redRef'][8]['MSE']
        # !! several runs; should be gathered together

        strings = self.session.allRes[0]['t1opt']['stringsTol']
        scores = self.session.allRes[0]['t1opt']['stringsTolScores']

        # ! reactions here is different. it should include the
        # AND edges as well
        if strings.ndim == 1:
            # BUGGY CNORfuzzy looks like bitstrings do not have correct length
            # if not enough strings are found.
            bstring = self.session.allRes[0]['t1opt']['bString']
            reactions = self.session.allRes[0]['processedModel']['reacID']
            N = len(bstring)
            M = len(reactions)
            fuzreactions = ['a=' + str(i) for i in range(0, N)]
            for i, reac in enumerate(reactions):
                fuzreactions[i] = reac
            df = pd.DataFrame([[0]*N], columns=fuzreactions)
        else:
            # FIXME what is it ? why adding a= ? reactions
            fuzreactions = ['a=' + str(i) for i in range(0, len(strings[0]))]
            for i, reac in enumerate(reactions):
                fuzreactions[i] = reac
            df = pd.DataFrame(strings, columns=fuzreactions)

        models = FuzzyModels(df)
        models.scores = scores
        models.cnograph.midas = self.data.copy()

        self.species = species
        self.reactions = reactions
        self.signals = self.session.signals

        # transforms the results into dataframes
        for i, res in enumerate(allRes):
            df = pd.DataFrame(res['t1opt']['results'],
                columns=("Generation","Best_score","Best_bitString","Stall_Generation",
                "Avg_Score_Gen","Best_score_Gen","Best_bit_Gen","Iter_time"))
            allRes[i]['t1opt']['results'] = df

        results = {
            'best_score': self.session.best_score,
            'best_bitstring': self.session.best_bitstring,
            'species': species,
        }

        self.results = FuzzyResults()
        self.results.results = results
        self.results.models = models
        self.results.allRes = allRes

    def _compute_mean_mses(self):
        """plot MSEs using interpolation of the results provided by the Fuzzy Analysis"""

        allFinalNumParams = self.session.summary['allFinalNumParams']
        dimRow = allFinalNumParams.shape[0]

        thresholds = self.thresholds[:]
        N = len(thresholds)
        allFinalMSEs = self.session.summary['allFinalMSEs']

        catExplore = np.zeros((dimRow, N))
        AbsNumParams = np.zeros((dimRow, N))
        AbsMSEs = np.zeros((dimRow, N))

        # interpolation
        for i in range(0, N):
            for j in range(0, dimRow):
                currIX = np.where(allFinalMSEs[j,] - allFinalMSEs[j,1] <= thresholds[i])[0]
                catExplore[j,i] = np.max(currIX)

        self.catExplore = catExplore

        for i in range(0,dimRow):
            for j in range(0, N):
                AbsNumParams[i,j] = allFinalNumParams[i, catExplore[i,j]]
                AbsMSEs[i,j] = allFinalMSEs[i, catExplore[i,j]]

        self.AbsMSEs = AbsMSEs
        self.AbsNumParams = AbsNumParams

        # final mean MSEs and number of parameters
        self.meanMSEs = np.mean(AbsMSEs, axis=0)
        self.meanNPs = np.mean(AbsNumParams, axis=0)

    def plot_mses(self, fontsize=20, **kwargs):
        """

        .. todo:: fix the yaxis and legend
        """
        self._compute_mean_mses()

        fig1 = pylab.figure()
        ax1 = fig1.add_subplot(111)

        line1 = ax1.semilogx(self.thresholds, self.meanMSEs, 'b-o', **kwargs)

        pylab.ylabel("MSEs", fontsize=fontsize)
        pylab.xticks(fontsize=16)
        pylab.yticks(fontsize=16)

        pylab.axis([self.thresholds[0], self.thresholds[-1],
        min(self.meanMSEs),max(self.meanMSEs)])

        ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)
        line2 = ax2.plot(self.thresholds, self.meanNPs, 'r-o')
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        pylab.ylabel("Number of Parameters", fontsize=fontsize)

        pylab.plot([], 'b-o', label="mean MSEs",)
        pylab.plot([], 'r-o', label="mean Number of Parameters")
        pylab.legend(loc="best", prop={'size':13})
        pylab.grid()
        pylab.xticks(fontsize=16)
        pylab.yticks(fontsize=16)

    def create_report_images(self):

        model = self.cnograph.copy()

        model.plot(filename=self._report._make_filename("pknmodel.svg"), show=False)
        model.preprocessing()
        model.plot(filename=self._report._make_filename("expmodel.png"), show=False)
        #self.plot_optimised_model(filename=self._make_filename("optimised_model.png"),
        #                          show=False)
        self.plot_errors(show=False)
        self._report.savefig("Errors.png")
        pylab.close()

        self.plot_mses(0.01)
        self._report.savefig("mse_vs_size.png")
        pylab.close()

        self.midas.plot()
        self._report.savefig("midas.png")
        pylab.close()

        self.plot_fitness(show=False, save=True)

    def plot_fitness(self, show=True, save=False):
        df = pd.DataFrame([res['t1opt']['results']['Best_score'].values
            for res in self.results.allRes])

        df = df.astype(float)

        pylab.clf()
        for res in self.results.allRes:
            pylab.plot(res['t1opt']['results']['Best_score'], '--', color='grey')
        pylab.grid()
        pylab.xlabel("Generation")
        pylab.ylabel("Score")
        #pylab.plot(df.mean().values, 'kx--', lw=3, label='Mean Score')

        y = df.mean().values
        x = range(0, len(y))
        yerr = df.std().values
        pylab.errorbar(x, y, yerr=yerr, xerr=None, fmt='-', label='Mean Score',
                color='k', lw=3)
        pylab.legend()

        if save is True:
            self._report.savefig("fitness.png")

        if show is False:
            pylab.close()

    def plot_errors(self, threshold=0.01, show=False):
        """

        :param threshold: fixed to 0.01 but should be provided by looking at :meth:`plot_mses`
        :param show:
        :return:
        """
        # TODO: should use automatic selection of the threshold
        # First, we need to compute the simulation
        script = """simulation = plotMeanFuzzyFit(%(threshold)s,
        summary$allFinalMSEs, allRes)""" % {'threshold': threshold}
        self.session.run(script)
        self.simulated = self.session.simulation['simResults']

        midas = self.data.copy()

        t0 = self.simulated['t0']
        t0 = pd.DataFrame(t0, columns=self.signals)
        t0.sort_index(axis=1, inplace=True)
        t0['experiment'] = midas.experiments.index
        t0['time'] = midas.times[0]
        t0['cell'] = midas.cellLines[0]

        t1 = self.simulated['t1']
        t1 = pd.DataFrame(t1, columns=self.signals)
        t1.sort_index(axis=1, inplace=True)
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
        try:
            midas.plot(mode="mse")
        except:
            pass
        return midas

    def create_report(self):

        self._create_report_header()

        txt = """<a href="http://www.cellnopt.org/">
            <object data="pknmodel.svg" type="image/svg+xml">
            <span>Your browser doesn't support SVG images</span> </object></a>"""
        txt += """<a class="reference external image-reference" href="scripts/exercice_3.py">
<img alt="MIDAS" class="align-right" src="midas.png" /></a>"""

        self._report.add_section(txt, "PKN graph", [("http://www.cellnopt.org", "cnograph")])

        self._report.add_section('<img src="expmodel.png">', "Expanded before optimisation")
        self._report.add_section( """<img src="optimised_model.png">""", "Optimised model")

        self._report.add_section('<img src="mse_vs_size.png">', "MSE vs Size")
        self._report.add_section('<img src="fitness.png">', "Fitness")
        self._report.add_section('<img src="Errors.png">', "Errors")

        #self._report.add_section(self.get_html_reproduce(), "Reproducibility")
        fh = open(self._report.directory + os.sep + "rerun.py", 'w')
        fh.write("from cno import CNORfuzzy\n")
        fh.write("CNORfuzzy(config=config.ini)\n")
        fh.write("c.optimise()\n")
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
        self._report.write("index.html")

    def _get_stats(self):
        res = {}
        #res['Computation time'] = self.total_time
        try:
            res['Best Score'] = self.results.results.best_score
        except:
            pass
        return res

    def _check_parameters(self, bs):
        if 'blength' not in self.__dict__:
            self._update()
        assert len(bs) == self.blength

    def _update(self):

        script = """
        library(CNORfuzzy)
        model = preprocessing(cnolist, pknmodel)
        indexlist = indexFinder(cnolist, model)
        params_fuzzy = defaultParametersFuzzy(cnolist, model)
        # model and cnolist reset internally if params provided
        simlist_fuzzy = prep4simFuzzy(model, params_fuzzy)"""
        self.session.run(script)

        simlist = self.session.simlist_fuzzy
        self.numType2 = simlist['numType2']
        self.numType1 = simlist['numType1']

        # numType1 are the number of edges from stimuli
        # numtype2 are the other edges
        # includes the AND gate edges.
        self.blength = simlist['numType2'] + simlist['numType1']

        params = self.session.params_fuzzy
        self.nTF = params['type1Funs'].shape[0]  # ==>seems to be type2Funs in computeScoreT1

    def create_random_parameters(self, update=True):
        import random
        if update is True:
            self._update()
        bs = [random.randint(0,self.nTF+1) for x in range(0,self.blength)]
        return bs

        #  intString <- (sample.int(dim(paramsList$type2Funs)[1],
        #         (simList$numType1+simList$numType2),replace=TRUE)) - 1
        #i d'abord numtype1

    def simulate(self, bstring, NAFac=1, sizeFac=0.0001):
        self._check_parameters(bstring)
        self.session.bstring = bstring
        params = {'NAFac': NAFac, 'sizeFac': sizeFac}
        script = """
        model = preprocessing(cnolist, pknmodel)
        indexlist = indexFinder(cnolist, model)
        params_fuzzy = defaultParametersFuzzy(cnolist, model)
        # model and cnolist reset internally if params provided
        simlist_fuzzy = prep4simFuzzy(model, params_fuzzy)
        score = computeScoreFuzzy(cnolist, model, simlist_fuzzy,
                indexlist, params_fuzzy,
                bstring, sizeFac=%(sizeFac)s, NAFac=%(NAFac)s)
        """ % params
        self.session.run(script)
        return self.session.score

def standalone(args=None):
    """This function is used by the standalone application called cellnopt_boolean

    ::

        cellnopt_boolean --help

    """
    if args is None:
        args = sys.argv[:]

    from cno.core.standalone import Standalone
    user_options = OptionsFuzzy()
    stander = Standalone(args, user_options)

    # just an alias
    options = stander.options

    if options.onweb is True or options.report is True:
        trainer = CNORfuzzy(options.pknmodel, options.data, verbose=options.verbose,
            verboseR=options.verboseR, config=options.config_file,  
            use_cnodata=options.cnodata)
        trainer.preprocessing()
    else:
        stander.help()

    params = stander.user_options.config.GA.as_dict()
    params.update(stander.user_options.config.Fuzzy.as_dict())
    trainer.optimise(**params)

    stander.trainer = trainer
    stander.report()


class OptionsFuzzy(OptionsBase):
    def __init__(self):
        prog = "cno_fuzzy"
        version = prog + " v1.0 (Thomas Cokelaer @2014)"
        super(OptionsFuzzy, self).__init__(version=version, prog=prog)
        self.add_section(ParamsGA())
        self.add_section(ParamsFuzzy())

if __name__ == "__main__":
    standalone()
