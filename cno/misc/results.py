from cno.misc.models import Models

import pylab
import pandas as pd


class Results(object):
    def __init__(self, name='default'):
        self.models = Models()

    def _get_scores(self):
        return self.models['score']
    scores = property(_get_scores)

    def _get_scores_vs_model_size_df(self):
        df = pd.DataFrame()
        df['mse'] = self.scores
        df['size'] = self.size
        return df

    def plot_fit(self):
        self.results.cnorbool.results[['Best_score','Avg_Score_Gen']].plot()
        pylab.title("Score per generation")
        pylab.xlabel("Generation")
        pylab.ylabel("Score")

    def hist_mse(self, fontsize=16, **kargs):
        """Plot histogram of the MSEs

         .. plot::
             :include-source:

             >>> from cellnopt.optimiser import ASPBool
             >>> from cellnopt.data import cnodata
             >>> a = ASPBool(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
             >>> a.run(fit=1)
             >>> a.hist_mse()


        """
        pylab.clf()
        mses = self.scores
        opt = self.scores.min()
        N = len(set(mses))
        print("There are %s different MSE found amongst %s models" % (N,len(mses)))
        res = pylab.hist(mses, **kargs)
        pylab.title("MSEs Distribution of the %s best models " % len(mses),
                    fontsize=fontsize)
        pylab.grid()
        pylab.plot([opt,opt], [0,max(res[0])], "r--",lw=2)
        pylab.xlabel("Mean Square Error of all models", fontsize=fontsize)
        pylab.ylabel("#",  fontsize=fontsize)

    def hist2d_scores_vs_model_size(self, bins=None, cmap='gist_heat_r',fontsize=16):
        from biokit import viz
        df = self._get_scores_vs_model_size_df()
        h = viz.Hist2d(df)
        if bins == None:
            mse_range = 20
            size_range = df.size.max() - df.size.min() + 1
            bins = (mse_range, size_range)
        h.plot(bins=bins, Nlevels=10, cmap=cmap)
        pylab.xlabel("Score", fontsize=fontsize)
        pylab.ylabel("Model size", fontsize=fontsize)

    def scatter_scores_vs_model_size(self):
        from biokit import viz
        df = self._get_scores_vs_model_size_df()
        viz.scatter_hist(df)


class BooleanResults(Results):
    def __init__(self):
        pass


class FuzzyResults(Results):
    def __init__(self):
        pass


class ODEResults(Results):
    def __init__(self):
        pass

