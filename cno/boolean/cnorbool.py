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

from easydev import Logging, AttrDict

from cno.io.multigraph import CNOGraphMultiEdges
from cno import CNOGraph, XMIDAS
from cno.core import CNOBase
from cno.misc.results import BooleanResults

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
    def __init__(self, model, data, verbose=True, verboseR=False):
        super(CNORbool, self).__init__(model, data, verbose=verbose)

        # this code could be moved to CNOBase
        self._verboseR = verboseR

        self.session = RSession(verbose=self.verboseR)
        self.parameters = {} # fill with GA binary parameters

    def _get_verboseR(self):
        return self._verboseR
    def _set_verboseR(self, value):
        self._verboseR = value
        self.session.dump_stdout = value
    verboseR = property(_get_verboseR, _set_verboseR)

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

        from cno.misc.models import Models
        df = pd.DataFrame(self.session.all_bitstrings,
                              columns=self.session.reactions)

        # cnograph created automatically from the reactions
        models = Models(df)
        models.cnograph.midas = self.data.copy()
        models.scores = self.session.all_scores

        results = {
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
        results['pkn'] = self.pknmodel
        results['midas'] = self.data
        self.results = BooleanResults()
        self.results.add_results(results)
        self.results.add_models(models)


    def plot_errors(self, close=False):

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
        if bs == None:
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


