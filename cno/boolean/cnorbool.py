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
import tempfile
import subprocess

from easydev import Logging, AttrDict

from cno.io.multigraph import CNOGraphMultiEdges
from cno import CNOGraph, XMIDAS
from cno.core import CNOBase

import pandas as pd
import pylab

__all__ = ["CNORbool"]

from biokit.rtools.session import RSession
from biokit.rtools import bool2R


class CNORbool(CNOBase):
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
    def __init__(self, model, data, verbose=True):
        super(CNORbool, self).__init__(model, data, verbose)
        self.verboseR = verbose
        self.session = RSession(dump_stdout=self.verboseR)
        self.parameters = {} # fill with GA binary parameters

    def reconnect(self):
        """If you cancel a job, you may need to reconnect the R server.

        """
        self.session = RSession(dump_stdout=self.verboseR)

    def optimise(self, tag="cnorbool", reltol=0.1, 
            expansion=True, compression=True):

        script_template = """
        library(CellNOptR)
        pknmodel = readSIF("%(pkn)s")
        cnolist = CNOlist("%(midas)s")
        #cnolist = normaliseCNOlist(cnolist, mode="ctrl", changeTh=40)
        model = preprocessing(cnolist, pknmodel, compression=%(compression)s, 
            expansion=%(expansion)s, maxInputsPerGate=3)

        res = gaBinaryT1(cnolist, model, relTol=%(reltol)s)
        sim_results = cutAndPlot(cnolist, model, list(res$bString),
                                 plotParams=list(maxrow = 80, cex=0.5),
                                 plotPDF=T, tag="%(tag)s")

        signals = colnames(cnolist@signals$`0`)
        colnames(sim_results$mse) = signals
        for (i in seq_along(sim_results$simResults[[1]])){
            colnames(sim_results$simResults[[1]][[i]]) = signals
        }

        #plotModel(model, cnolist, bString=res$bString,
        #          filename="Model-CNO-%(tag)s", output="PDF")

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
        self.results = {}
        self.results['cnorbool'] = {}

        script = script_template % {
                'pkn': self.pknmodel.filename, 
                'midas': self.data.filename,
                'tag':tag,
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

        from cno.misc.models import Models
        df = pd.DataFrame(self.session.all_bitstrings,
                              columns=self.session.reactions)
        models = Models(df)
        models.cnograph.midas = self.data.copy()

        self.results['cnorbool'] = {
                'best_score': self.session.best_score,
                'best_bitstring': self.session.best_bitstring,
                'all_scores': self.session.all_scores,
                'all_bitstrings': self.session.all_bitstrings,
                'reactions': self.session.reactions,
                'sim_results': self.session.sim_results,  # contains mse and sim at t0,t1,
                'reactions': self.session.reactions,
                'results': results,
                'models':models,
                'stimuli':self.session.stimuli.copy(),
                'inhibitors':self.session.inhibitors.copy(),
                'species':self.session.species,
                'tag': tag
        }
        self.results['pkn'] = self.pknmodel
        self.results['midas'] = self.data
        self.results = AttrDict(**self.results)
        self.results['cnorbool'] = AttrDict(**self.results['cnorbool'])

    def plot_errors(self, close=False):

        results = self.results['cnorbool']
        midas = self.results['midas']

        t0 = results['sim_results']['simResults'][0]['t0']
        t1 = results['sim_results']['simResults'][0]['t1']
        mse = results['sim_results']['mse']
        stimuli = results['stimuli']
        inhibitors = results['inhibitors']
        species = results['species']
        tag = results['tag']

        # need to make sure that order in cellnopt is the same as in the midas
        # shoudl be fine
        assert all([all(midas.stimuli.ix[0] == stimuli.ix[0]) for i in range(0,40)])
        assert all([all(midas.inhibitors.ix[0] == inhibitors.ix[0]) for i in range(0,40)])

        t0 = pd.DataFrame(t0, columns=species)
        t0['experiment'] = midas.experiments.index
        t0['time'] = midas.times[0]
        t0['cellLine'] = midas.cellLines[0]

        t1 = pd.DataFrame(t1, columns=species)
        t1['experiment'] = midas.experiments.index
        t1['time'] = midas.times[1]
        t1['cellLine'] = midas.cellLines[0]

        df = pd.concat([t0,t1]).set_index(['cellLine', 'experiment', 'time'])
        df.sortlevel(1, inplace=True)

        midas.sim = df.copy()
        midas.cmap_scale = 1   # same a CellNOptR
        midas.plot(mode="mse")

        #midas.plotSim()
        pylab.savefig("Error-{0}.png".format(tag), dpi=200)
        pylab.savefig("Error-{0}.svg".format(tag), dpi=200)
        if close:
            pylab.close()
        return midas

    def __str__(self):
        txt = ""
        for cell in self.cellLines:
            txt += "{0}: {1}".format(cell, self.results['best_score']) + "\n"
        return txt

    def simulate(self):
        """


        c.results.cnorbool.best_bitstring
        array([1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0])
        c.results.cnorbool.reactions


        """


        script = """

        library(CellNOptR)

        """

        self.session.run(script)
       
    def _get_models(self):
        return self.results.cnorbool.models
    models = property(_get_models)
