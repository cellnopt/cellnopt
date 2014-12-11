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

from easydev import Logging, AttrDict

from cno.io.multigraph import CNOGraphMultiEdges
from cno.core import CNOBase, CNORBase
from cno.core.results import BooleanResults
from cno.core.models import BooleanModels
from cno.core import ReportBool

import pandas as pd
import pylab
from biokit.rtools import bool2R


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
    def __init__(self, model, data, tag=None, verbose=True, verboseR=False,
            config=None):
        CNOBase.__init__(self,model, data, tag=tag, verbose=verbose,
                config=config)
        CNORBase.__init__(self, verboseR)
        self.parameters = {} # fill with GA binary parameters
        self._report = ReportBool()
        self.results = BooleanResults()

    def optimise(self, tag="cnorbool", reltol=0.1,
            expansion=True, maxgens=150, stallgenmax=100, compression=True):

        script_template = """
        library(CellNOptR)
        pknmodel = readSIF("%(pkn)s")
        cnolist = CNOlist("%(midas)s")
        model = preprocessing(cnolist, pknmodel, compression=%(compression)s,
            expansion=%(expansion)s, maxInputsPerGate=3)

        res = gaBinaryT1(cnolist, model, relTol=%(reltol)s,
            stallGenMax=%(stallgenmax)s, maxGens=%(maxgens)s)
        sim_results = cutAndPlot(cnolist, model, list(res$bString),
                                 plotParams=list(maxrow = 80, cex=0.5),
                                 plotPDF=T, tag="%(tag)s")

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
        self._results = {}
        script = script_template % {
                'pkn': self.pknmodel.filename,
                'midas': self.data.filename,
                'tag':tag,
                'maxgens':maxgens,
                'stallgenmax':stallgenmax,
                'reltol':reltol,
                'compression': bool2R(compression),
                'expansion': bool2R(expansion)
                }
        self.session.run(script)

        # need to change type of some columns, which are all string
        results = self.session.results

        columns_int = ['Generation', 'Stall_Generation']
        columns_float = ['Best_score', 'Avg_Score_Gen', 'Best_score_Gen', 'Iter_time']
        results[columns_int] = results[columns_int].astype(int)
        results[columns_float] = results[columns_float].astype(float)

        # cnograph created automatically from the reactions
        df = pd.DataFrame(self.session.all_bitstrings,
                              columns=self.session.reactions)
        models = BooleanModels(df)
        models.scores = self.session.all_scores
        models.cnograph.midas = self.data.copy()

        reactions = self._reac_cnor2cno(self.session.reactions)

        results = {
                'best_score': self.session.best_score,
                'best_bitstring': self.session.best_bitstring,
                'all_scores': self.session.all_scores,
                'all_bitstrings': self.session.all_bitstrings,
                'reactions': reactions,
                'sim_results': self.session.sim_results,  # contains mse and sim at t0,t1,
                'results': results,
                'models':models,
                'stimuli':self.session.stimuli.copy(),
                'inhibitors':self.session.inhibitors.copy(),
                'species':self.session.species,
                'tag': tag
        }
        results['pkn'] = self.pknmodel
        results['midas'] = self.data

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
        for node in self._pknmodel.nodes():
            self._pknmodel.node[node]['URL'] = "http://www.uniprot.org/uniprot/?query=Ras&sort=score"


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

    def onweb(self, show=True):
        self.create_report()
        try:
            self.create_report_images()
        except:
            pass
        if show:        self._report.show()

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
        txt += "o.report()\n</pre>\n"
        self._report.add_section(txt, "Script used")

        txt = """<a href="http://www.cellnopt.org/">
            <object data="pknmodel.svg" type="image/svg+xml" >
            <span>Your browser doesn't support SVG images</span>

            </object></a>"""

        # <img src="pknmodel.png">

        self._report.add_section(txt, "PKN graph", [("http://www.cellnopt.org", "cnograph")])

        txt = """<a class="reference external image-reference" href="scripts/exercice_3.py">
<img alt="MIDAS" class="align-left" src="midas.png" /></a>"""
        self._report.add_section(txt, "Data", [("http://www.cellnopt.org", "XMIDAS class")])

        self._report.add_section('<img class="figure" src="expmodel.png">',
            "Expanded before optimisation")
        self._report.add_section('img class="figure" src="optimised_model.png">',
            "Optimised model")
        self._report.add_section('<img class="figure" src="optimised_model_mapback.png">',
            "Optimised mapback model")
        self._report.add_section('<img class="figure" src="fitness.png">',
            "Fitness")
        self._report.add_section('<img class="figure" src="Errors.png">',
            "Errors")

        self._report.add_section(self.get_html_reproduce(), "Reproducibility")
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
            txt += "<tr><td>%s</td><td>%s</td></tr>\n" % (k, stats[k])
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

