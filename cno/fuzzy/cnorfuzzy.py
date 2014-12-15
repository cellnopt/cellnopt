import sys
import os

import pandas as pd
import numpy as np
import pylab

from cno.core import CNOBase
from cno.core.params import OptionsBase

from cno.core import CNORBase
from cno.misc.results import FuzzyResults
from cno.core.report import ReportFuzzy
from cno.core.params import ParamsGA
from biokit.rtools import bool2R


__all__ = ["CNORfuzzy"]

"""
$redThresh  0e+00 1e-04 5e-04 1e-03 3e-03 5e-03 1e-02
$doRefinement  TRUE
$sizeFac  1e-04
$NAFac  1
$popSize  50
$pMutation  0.5
$maxTime  180
$maxGens  500
$stallGenMax  100
$selPress  1.2
$elitism  5
$relTol  0.1
$verbose  TRUE
$optimisation
    $optimisation$algorithm "NLOPT_LN_SBPLX"
    $optimisation$xtol_abs  0.001
    $optimisation$maxEval  10000
    $optimisation$maxTime  300

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
    def __init__(self, model=None, data=None, verbose=True, verboseR=False):
        """.. rubric:: constructor"""
        #super(CNORfuzzy, self).__init__(model, data, verbose=verbose)
        CNOBase.__init__(self,model, data, verbose=verbose)
        CNORBase.__init__(self, verboseR=verboseR)

        self.parameters = FuzzyParameters()

        self._report = ReportFuzzy()
        self._report._init_report()

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

    def optimise(self, tag="cnorfuzzy", N=2,
            popsize=50,reltol=0.1, maxtime=180, expansion=True, maxgens=150,
            stallgenmax=100, compression=True):


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
        #for x in mapping.keys():
        #    # update config data structure only if user parameter provided
        #    if x in kargs.keys():
        #        self.config.GA[mapping[x]] = kargs[x]
        #    params[x] = self.config.GA[mapping[x]]
        #
        #    elitism=%(elitism)s, pMutation=%(pmutation)s,
        #    NAFac=%(nafac)s,  selPress=%(selpress)s, relTol=%(reltol)s, sizeFac=%(sizefac)s,
        #    stallGenMax=%(stallgenmax)s)
        script_template = """
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
        paramsList$stallGenMax = %(stallgenmax)s
        paramsList$optimisation$maxtime = 60*5

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
        """
        script = script_template % {
                'pkn': self.pknmodel.filename,
                'midas': self.data.filename,
                'tag': tag,
                'N': N,
                'popsize': popsize,
                'maxgens': maxgens,
                'maxtime': maxtime,
                'stallgenmax': stallgenmax,
                'reltol': reltol,
                'compression': bool2R(compression),
                'expansion': bool2R(expansion)
                }
        self.session.run(script)
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

        #scores = res1['t1opt']['stringsTolScores']
        # res1['t1opt']['stringsTol']

        reactions = res1['processedModel']['reacID']
        bScore = res1['currBestDiscrete']
        bString = res1['intString']

        # TODO: find best MSE/bitstring amongst the N runs

        # redRef contains a MSE for each threshold
        res1['redRef'][8]['MSE']

        # !! several runs; should be gathered together
        from cno.misc.models import FuzzyModels

        strings = self.session.allRes[0]['t1opt']['stringsTol']
        scores = self.session.allRes[0]['t1opt']['stringsTolScores']

        # ! reactions here is different. it should include the
        # AND edges as well
        print(len(reactions), len(strings[0]))
        fuzreactions = ['a=' + str(i) for i in range(0, len(strings[0]))]
        for i, reac in enumerate(reactions):
            fuzreactions[i] = reac
        df = pd.DataFrame(strings, columns=fuzreactions)
        models = FuzzyModels(df)
        models.scores = scores
        models.cnograph.midas = self.data.copy()

        self.species = species
        self.reactions = reactions
        self.bScore = bScore
        self.bString = bString
        self.signals = self.session.signals

        # transforms the results into dataframes
        for i, res in enumerate(allRes):
            df = pd.DataFrame(res['t1opt']['results'],
                columns=("Generation","Best_score","Best_bitString","Stall_Generation",
                "Avg_Score_Gen","Best_score_Gen","Best_bit_Gen","Iter_time"))
            allRes[i]['t1opt']['results'] = df

        results = {}
        results['species'] = species
        results['bScore'] = bScore

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

    def todo(self):
        from cno.core.models import Models
        res = self.results.allResp[0]
        Models( pd.DataFrame(res.t1opt['stringsTol'], columns=list(self.reactions)))

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
        txt += "todo\n"
        txt += "o.report()\n</pre>\n"
        self._report.add_section(txt, "Script used")

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

        self._report.add_section(self.get_html_reproduce(), "Reproducibility")
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
        self._report.write(self._report.directory, "index.html")

    def _get_stats(self):
        res = {}
        #res['Computation time'] = self.total_time
        try:
            res['Best Score'] = self.bScore
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


        #
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

        cno_fuzzy --help

    """
    if args == None:
        args = sys.argv[:]

    user_options = OptionFuzzy(prog="cno_fuzzy")

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


class OptionFuzzy(OptionsBase):

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



