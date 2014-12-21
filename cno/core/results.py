
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
from cno.core.models import Models

from biokit import viz

import pylab
import pandas as pd

from easydev import Logging, AttrDict


__all__ = ['BooleanResults', 'FuzzyResults', 'ODEResults', 'DTResults']


class Results(object):
    """

    Plots whatever is related to the size of the models (e.g., number of active edges)
    and the score (e.g., MSE)

    These are the most generic information thta could be retrieved and that are not specific
    to a formalism.

    For instance, histogram of the scores.
    """
    def __init__(self):
        self._models = None
        self._results = None

    def _get_models(self):
        return self._models
    def _set_models(self, models):
        self._models = models.copy()
    models = property(_get_models, _set_models)

    def _get_results(self):
        return self._results
    def _set_results(self, results):
        self._results = AttrDict(**results)
    results = property(_get_results, _set_results)

    def _get_scores(self):
        return self.models.scores
    scores = property(_get_scores)

    def _get_sizes(self):
        return self.models.sizes
    sizes = property(_get_sizes)


class BooleanResults(Results):

    def __init__(self):
        pass

    def plot_fit(self):
        self.results.results[['Best_score','Avg_Score_Gen']].plot()
        pylab.title("Score per generation")
        pylab.xlabel("Generation")
        pylab.ylabel("Score")

    def _get_scores_vs_model_size_df(self):
        df = pd.DataFrame()
        df['score'] = self.scores
        df['size'] = self.sizes
        return df

    def hist2d_scores_vs_model_size(self, bins=None, cmap='gist_heat_r', 
            fontsize=16, contour=False, Nlevels=10):
        df = self._get_scores_vs_model_size_df()
        h = viz.Hist2d(df)
        if bins == None:
            scores_range = 20
            size_range = df.size.max() - df.size.min() + 1
            bins = (scores_range, size_range)
        h.plot(bins=bins, Nlevels=Nlevels, cmap=cmap, contour=contour)
        pylab.xlabel("Score", fontsize=fontsize)
        pylab.ylabel("Model size", fontsize=fontsize)

    def scatter_scores_vs_model_size(self):
        """Scatter plot of the model size and scores"""
        df = self._get_scores_vs_model_size_df()
        viz.scatter_hist(df)

    def hist_scores(self, fontsize=16, **kargs):
        """Plot histogram of the MSEs

         .. plot::
             :include-source:

             >>> from cno.boolean.cnorbool import CNORbool
             >>> from cno import cnodata
             >>> a = CNORbool(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
             >>> a.optimise()
             >>> a.results.hist_scores()

        """
        pylab.clf()
        scores = self.scores
        opt = self.scores.min()
        N = len(set(scores))
        print("There are %s different MSE found amongst %s models" % (N,len(scores)))
        res = pylab.hist(scores, **kargs)
        pylab.title("MSEs Distribution of the %s best models " % len(scores),
                    fontsize=fontsize)
        pylab.grid()
        pylab.plot([opt,opt], [0,max(res[0])], "r--",lw=2)
        pylab.xlabel("Mean Square Error of all models", fontsize=fontsize)
        pylab.ylabel("#",  fontsize=fontsize)


class FuzzyResults(Results):
    def __init__(self):
        pass


class DTResults(BooleanResults):
    def __init__(self):
        super(DTResults, self).__init__()


class ODEResults(Results):
    def __init__(self):
        pass

    def plot_fit(self, ymin=None, ymax=None):
        pylab.plot(self.results.all_scores, 'o-')
        pylab.grid()
        pylab.title("Best score per generation")
        pylab.xlabel("Generation")
        pylab.ylabel("Score")
        ylim = pylab.ylim()
        if ymin is None:
            ymin = 0
        if ymax is None:
            ymax = ylim[1] * 1.5
        pylab.ylim([ymin, ymax])
        pylab.grid(True)

