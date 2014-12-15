# -*- python -*-
#
#  This file is part of CNO software
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
import sys

import numpy as np
import pandas as pd
import pylab

from cno.core import CNOBase, CNORBase
from cno.core.results import BooleanResults
from cno.core.models import BooleanModels
from cno.core import ReportBool
from cno.core.params import OptionsBase

from biokit.rtools import bool2R

from cno.core.params import params_to_update



__all__ = ["CNORbool"]


class CNORbool(CNOBase, CNORBase):
    """Access to CellNOptR R package to run boolean analysis


    ::

        c = pipeline.CNObool("PKN-test.sif", "MD-test.csv")
        c.optimise(compression=True, expansion=True, reltol=.15)


    Results are stored in :attr:`results`. Information stored are various.
    The errors corresponding to the best models can be visualised with :meth:`plot_errors`
    and models within the tolerance are stored in :attr:`models.

    .. plot::
        :include-source:

        from cno import cnodata, CNORbool
        c = CNORbool(cnodata("PKN-ToyMMB.sif"),
            cnodata("MD-ToyMMB.csv"))
        c.optimise()
        c.plot_errors()

    """
    #params = BooleanParameters.default
    def __init__(self, model, data, tag=None, verbose=True,
                 verboseR=False, config=None):

        CNOBase.__init__(self,model, data, tag=tag, verbose=verbose,
                config=config)
        CNORBase.__init__(self, verboseR)

        self._report = ReportBool()
        self.results = BooleanResults()
        self.results2 = BooleanResults()

        self.config.General.pknmodel.value = self.pknmodel.filename
        self.config.General.data.value = self.data.filename

    def _update_config(self, section, kwargs):
        for k, v in kwargs.items():
            try:
                section = getattr(self.config, section)
                option = getattr(section, k)
                option.value = v
            except:
                # additional options should be ignored.
                pass

    # !! should be same default values as in paramsGA
    @params_to_update()
    def optimise(self,  NAFac=1, pmutation=0.5, selpress=1.2, popsize=50,
            reltol=0.1, elistim=5, maxtime=60, sizefactor=0.0001,
            time_index_1=1, maxgens=500, maxstallgens=100, another=True):

        # update config GA section with user parameters
        self._update_config('GA', self.optimise.actual_kwargs)

        # keep track of the GA parameters, which may have been update above
        gad = dict([(k, self.config.GA[k].value)
            for k in self.config.GA._get_names()])

        script_template = """
        library(CellNOptR)
        pknmodel = readSIF("%(pkn)s")
        cnolist = CNOlist("%(midas)s")
        model = preprocessing(cnolist, pknmodel, compression=%(compression)s,
            expansion=%(expansion)s, maxInputsPerGate=3)

        res = gaBinaryT1(cnolist, model, popSize=%(popsize)s, maxGens=%(maxgens)s,
            maxTime=%(maxtime)s, elitism=%(elitism)s, pMutation=%(pmutation)s,
            NAFac=%(NAFac)s,  selPress=%(selpress)s, relTol=%(reltol)s, sizeFac=%(sizefactor)s,
            stallGenMax=%(maxstallgens)s)

        sim_results = cutAndPlot(cnolist, model, list(res$bString),
                                 plotParams=list(maxrow = 80, cex=0.5),
                                 plotPDF=F)
        sim_results2 = NULL

        signals = colnames(cnolist@signals$`0`)
        colnames(sim_results$mse) = signals
        for (i in seq_along(sim_results$simResults[[1]])){
            colnames(sim_results$simResults[[1]][[i]]) = signals
        }

        # to be retrieved inside Python code
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

        params = {
            'pkn': self.pknmodel.filename,
            'midas': self.data.filename,
            'compression': bool2R(compression),
            'expansion': bool2R(expansion)
            }
        params.update(gad)

        script = script_template % params
        self.session.run(script)
        self.reactions_r = self.session.reactions
        # need to change type of some columns, which are all string
        results = self.session.results
        results.columns = [x.strip() for x in results.columns]

        columns_int = ['Generation', 'Stall_Generation']
        columns_float = ['Best_score', 'Avg_Score_Gen', 'Best_score_Gen', 'Iter_time']
        results[columns_int] = results[columns_int].astype(int)
        results[columns_float] = results[columns_float].astype(float)

        # cnograph created automatically from the reactions
        df = pd.DataFrame(self.session.all_bitstrings,
                              columns=self.reactions_r)
        models = BooleanModels(df)
        # flatten to handle exhaustive
        import numpy as np
        models.scores = np.array(list(pylab.flatten(self.session.all_scores)))
        models.cnograph.midas = self.data.copy()


        results = {
                'best_score': self.session.best_score,
                'best_bitstring': self.session.best_bitstring,
                'all_scores': self.session.all_scores,
                'all_bitstrings': self.session.all_bitstrings,
                'reactions': self._reac_cnor2cno(self.reactions_r),
                'sim_results': self.session.sim_results,  # contains mse and sim at t0,t1,
                'results': results,
                'models':models,
                'stimuli':self.session.stimuli.copy(),
                'inhibitors':self.session.inhibitors.copy(),
                'species':self.session.species,
        }
        results['pkn'] = self.pknmodel
        results['midas'] = self.data

        self.results.results = results
        self.results.models = models

    #use same parameters as in T1
    def optimise2(self, reltol=0.1, time_index_2=3, best_bitstring=None):
        # TODO assert there are 2 time indices
        # something interesting to do is to run steady stte not only
        # for the best bitstring but all those within the tolerance !

        if best_bitstring is None:
            best_bitstring = self.results.results.best_bitstring

        reactions = [r for r,b in zip(self.reactions_r, best_bitstring) if b==0]
        script_template = """
        # model, cnolist and best_bitstring are created when calling
        # optimise()
        res2 = gaBinaryTN(cnolist, model, popSize=%(popsize)s, maxGens=%(maxgens)s,
            maxTime=%(maxtime)s, elitism=%(elitism)s, pMutation=%(pmutation)s,
            NAFac=%(NAFac)s,  selPress=%(selpress)s, relTol=%(reltol)s, sizeFac=%(sizefactor)s,
            stallGenMax=%(maxstallgens)s, timeIndex=%(timeindex)s, bStrings=list(best_bitstring))
            # bitstring is provided when calling optimise()


        sim_results2 = cutAndPlot(cnolist, model, list(res$bString, res2$bString),
                                   plotParams=list(maxrow = 80, cex=0.5),
                                                                    plotPDF=F)


        results2 = as.data.frame(res2$results)
        best_bitstring2 = res2$bString
        best_score2 = res2$bScore
        all_scores2 = res2$stringsTolScores
        all_bitstrings2 = res2$stringsTol
        """
        params = dict([(k, self.config.GA[k].value)
            for k in self.config.GA._get_names()])
        params['timeindex'] = time_index_2
        params['reltol'] = reltol

        script= script_template % params
        self.session.run(script)

        results = self.session.results2
        results.columns = [x.strip() for x in results.columns]
        columns_int = ['Generation', 'Stall_Generation']
        columns_float = ['Best_score', 'Avg_Score_Gen', 'Best_score_Gen', 'Iter_time']
        results[columns_int] = results[columns_int].astype(int)
        results[columns_float] = results[columns_float].astype(float)

        # cnograph created automatically from the reactions
        # TODneed to gure outat are  reactio names

        try:
            print("only 1 bitstring")
            df = pd.DataFrame(self.session.all_bitstrings2,
                    columns=reactions)
            self.df = df
        except:
            print("several bitstrings")
            df = pd.DataFrame(self.session.all_bitstrings2)
            df = df.transpose()
            self.reactions2 = reactions
            df.columns = reactions
            self.df = df

        models = BooleanModels(df)
        # flatten to handle exhaustive
        try:
            # if only 1 item, it is not a list...
            models.scores = np.array(list(pylab.flatten(self.session.all_scores2)))
        except:
            models.scores = np.array([self.session.all_scores2])

        #models.cnograph.midas = self.data.copy()
        #reactions = self._reac_cnor2cno(self.session.reactions)

        results = {
                'best_score': self.session.best_score2,
                'best_bitstring': self.session.best_bitstring2,
                'all_scores': self.session.all_scores2,
                'all_bitstrings': self.session.all_bitstrings2,
                'reactions': reactions,
                'sim_results': self.session.sim_results,  # contains mse and sim at t0,t1,
                'results': results,
                'models':models,
                'stimuli':self.session.stimuli.copy(),
                'inhibitors':self.session.inhibitors.copy(),
                'species':self.session.species,
        }
        results['pkn'] = self.pknmodel
        results['midas'] = self.data

        self.results2.results = results
        self.results2.models = models


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

    def simulate(self, bs=None, compression=True, expansion=True):
        """

        input could be a bitstring with correct length and same order
        OR a model


        c.results.cnorbool.best_bitstring
        c.results.cnorbool.reactions
        array([1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0])
        c.results.cnorbool.reactions

        """
        if bs is None:
            bs = ",".join([str(x) for x in self.results.best_bitstring])
        else:
            bs = ",".join([str(x) for x in bs])

        script_template = """
        library(CellNOptR)
        pknmodel = readSIF("%(pkn)s")
        cnolist = CNOlist("%(midas)s")
        model = preprocessing(cnolist, pknmodel, compression=%(compression)s,
            expansion=%(expansion)s, maxInputsPerGate=3)
        mse = computeScoreT1(cnolist, model, %(bs)s)
        """

        script = script_template % {
                'pkn': self.pknmodel.filename,
                'midas': self.data.filename,
                'compression': bool2R(compression),
                'expansion': bool2R(expansion),
                'bs': "c(" + bs + ")"
                }

        self.session.run(script)
        return self.session.mse

    def _get_models(self):
        return self.results.models
    models = property(_get_models)

    def create_report_images(self):
        # ust a simple example of settinh the uniprot url
        # should be part of cellnopt.core
        #for node in self._pknmodel.nodes():
        #    self._pknmodel.node[node]['URL'] = "http://www.uniprot.org/uniprot/?query=Ras&sort=score"

        self._pknmodel.plot(filename=self._report._make_filename("pknmodel.svg"),
                            show=False)
        self._pknmodel.plot(filename=self._report._make_filename("pknmodel.png"),
                            show=False)

        model = self.cnograph.copy()
        model.preprocessing()
        model.plot(filename=self._report._make_filename("expmodel.png"), show=False)

        self.plot_optimised_model(filename=self._report._make_filename("optimised_model.png"),
                                  show=False)

        self.plot_mapback_model()
        self._report.savefig("optimised_model_mapback.png")

        self.plot_errors(show=False)
        self._report.savefig("Errors.png")

        self.midas.plot()
        self._report.savefig("midas.png")
        pylab.close()

        self.plot_fitness(show=False, save=True)
        # self._report.savefig("plot_fit.png")

    def plot_fitness(self, show=True, save=False):
        # TODO show and save parameters
        self.results.plot_fit()

        if save is True:
            self._report.savefig("fitness.png")

        if show is False:
            pylab.close()

    def create_report(self):
        self._report._init_report()
        self._report.directory = self._report.report_directory
        # Save filenames and report in a section
        fname = self._report.directory + os.sep + "PKN-pipeline.sif"

        self.config.save(self._report.directory + os.sep + 'config.ini')
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

        txt = """<a href="http://www.cellnopt.org/">
            <object data="pknmodel.svg" type="image/svg+xml" >
            <span>Your browser doesn't support SVG images</span>

            </object></a>"""

        # <img src="pknmodel.png">

        self._report.add_section(txt, "PKN graph", [("http://pythonhosted.org//cno/references.html#module-cno.io.cnograph", "cnograph")])

        txt = """<a class="reference external image-reference" href="scripts/exercice_3.py">
<img alt="MIDAS" class="align-left" src="midas.png" /></a>"""
        self._report.add_section(txt, "Data", [("http://pythonhosted.org//cno/references.html#module-cno.io.midas.xmidas", "XMIDAS class")])

        self._report.add_section('<img class="figure" src="expmodel.png">',
            "Expanded before optimisation")
        self._report.add_section('<img class="figure" src="optimised_model.png">',
            "Optimised model")
        self._report.add_section('<img class="figure" src="optimised_model_mapback.png">',
            "Optimised mapback model")
        self._report.add_section('<img class="figure" src="fitness.png">',
            "Fitness")
        self._report.add_section('<img class="figure" src="Errors.png">',
            "Errors")

        #self._report.add_section(self.get_html_reproduce(), "Reproducibility")
        fh = open(self._report.report_directory + os.sep + "rerun.py", 'w')
        fh.write("from cellnopt import CNORbool, cnodata*\n")
        fh.write("CNORbool(config=config.ini)\n")
        fh.write("c.optimise()\n")
        fh.write("c.onweb()\n")
        fh.close()

        # some stats
        stats = self._get_stats()
        txt = "<table>\n"
        for k in sorted(stats.keys()):
            txt += "<tr><td>%s:</td><td>%s</td></tr>\n" % (k, stats[k])
        txt += "</table>\n"
        #txt += """<img id="img" onclick='changeImage();' src="fit_over_time.png">\n"""
        self._report.add_section(txt, "stats")

        # dependencies
        self._report.add_section('<img src="dependencies.svg">', "Dependencies")
        self._report.write(self._report.report_directory, "index.html")

    def _get_stats(self):
        res = {}
        #res['Computation time'] = self.total_time
        try:
            res['Best Score'] = self.results.results.best_score
            res['Total time'] = self.results.results.results.Iter_time.sum()
        except:
            pass
        return res


def standalone(args=None):
    """This function is used by the standalone application called cellnopt_boolean

    ::

        cellnopt_boolean --help

    """
    if args == None:
        args = sys.argv[:]

    user_options = OptionsBoolean()

    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    else:
        options = user_options.parse_args(args[1:])

    if options.onweb is True or options.report is True:
        o = CNORbool(options.pknmodel, options.data, verbose=options.verbose,
                 verboseR=options.verboseR, config=user_options.config)

    if options.onweb is True:
        o.optimise()
        o.onweb()
    elif options.report is True:
        o.optimise()
        o.report()
    else:
        from easydev.console import red
        print(red("No report requested; nothing will be saved or shown"))
        print("use --on-web or --report options")

class OptionsBoolean(OptionsBase):

    def  __init__(self):
        prog = "cellnopt_boolean_steady"
        version = prog + " v1.0 (Thomas Cokelaer @2014)"
        super(OptionsBoolean, self).__init__(version=version, prog=prog)
        from cno.core.params import ParamsGA
        section = ParamsGA()
        self.add_section(section)





if __name__ == "__main__":
    """Used by setup.py as an entry point to :func:`standalone`"""
    standalone(sys.argv)



