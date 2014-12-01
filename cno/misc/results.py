from cno.misc.models import Models

from biokit import viz

import pylab
import pandas as pd


from easydev import Logging, AttrDict


class Results(object):
    """

    Plots whatever is related to the size of the models (e.g., number of active edges)
    and the score (e.g., MSE)

    These are the most generic information thta could be retrieved and that are not specific
    to a formalism.

    For instance, histogram of the scores.
    """
    def __init__(self, models):
        self.models = models.copy()

    def add_results(self, results):
         self.results = AttrDict(**results)
         #self.results[name] = AttrDict(**self.results[name])

    def add_models(self, models):
        pass

    def _get_scores(self):
        return self.models['score']
    scores = property(_get_scores)





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
        df['size'] = self.size
        return df

    def hist2d_scores_vs_model_size(self, bins=None, cmap='gist_heat_r', 
            fontsize=16):
        df = self._get_scores_vs_model_size_df()
        h = viz.Hist2d(df)
        if bins == None:
            scores_range = 20
            size_range = df.size.max() - df.size.min() + 1
            bins = (scores_range, size_range)
        h.plot(bins=bins, Nlevels=10, cmap=cmap)
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

             >>> from cellnopt.optimiser import ASPBool
             >>> from cellnopt.data import cnodata
             >>> a = ASPBool(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
             >>> a.run(fit=1)
             >>> a.hist_scores()

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



