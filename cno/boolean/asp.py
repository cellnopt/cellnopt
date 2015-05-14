import subprocess
import os
from collections import defaultdict
from collections import Counter

from cno.io import CNOGraph, XMIDAS
from cno.core.models import BooleanModels

from easydev import AttrDict, TempFile

import pylab
import pandas as pd
from biokit import viz

from cno.core.base import CNOBase


__all__ = ["CASPO"]


def combine_caspo(results):
    """From several results from CASPO, create a MultiCASPO object"""
    raise NotImplementedError


class MultiCASPO(object):
    def __init__(self, verbose=False):
        self.verbose = verbose


class CASPO(CNOBase):
    def __init__(self, pknmodel, data, verbose=False):
        super(CASPO, self).__init__(pknmodel, data, verbose)




        # TODO out and err should be temp file
        out = TempFile('out')
        err = TempFile('err')
        status = os.system('caspo --help 1>{0} 2>{1}'.format(out.name, err.name))
        if status !=0:
            raise ValueError("Could not find caspo executable")

        # TODO check if caspo executable is available and works!
        # TODO compression ??
        # TODO save results in a temporary directory

        self.df = pd.DataFrame({'mse':[], 'size':[]})

    def optimise(self, size=0, fit=0, factor=10, length=0, timepoint=None):
        """Run the optimisation using **clasp/gringo**.

        :param bool gtts: if set to True, in addition to the models, the GTTS are
             also stored. This could increase the computational time
            dramatically.
       :param float fit: this parameter allows user to obtain all models with a
             fit/MSE in the range [:math:`MSE_0`; :math:`MSE_0 \\times fit \\times
             100`] where :math:`MSE_0` is the best fit of the optimal model and
             :attr:`fit` is in percentage. Be aware that the number of models may
             increase dramatically as well as the computation time if the
             required fit is too large. THis is network/data dependent.
         :param int size: if the optimal model has size N, only the models with size
             smaller than N + size will be kept and returned.
         :param discrete: number of points used in the discretization. Could be
             used to decrease computation time


        """

        # save local version of the PKN and MIDAS (CASPO requires it...)
        from easydev import TempFile
        sifname = TempFile(suffix='.sif')
        self.pknmodel.to_sif(sifname.name)
        midasname = TempFile(suffix='.csv')
        self.midas.to_midas(expand_time_column=True, filename=midasname.name)

        if timepoint == None:
            timepoint = self.midas.times[1]

        # TODO: discretisation can be floor, ceil, floor
        # TODO length option (max length for conjunecitons default to 0 ?)

        # FIXME: bug in caspo. when fit=0, actually it is set to size of models +fit + 1
        #   so internally, we set fit to fit-1
        # size -= 1

        script = "caspo learn %(pkn)s %(midas)s %(timepoint)s --fit %(fit)s --size %(size)s --factor %(factor)s --discretization %(disc)s --length %(length)s" % \
                 {'midas': midasname.name, 'pkn': sifname.name, 'fit':fit, 'length':length,
        'size':size, 'factor':factor, 'timepoint':timepoint, 'disc':'round'}
        self.logging.info(script)

        ret = subprocess.Popen("%s" % (script),
            stdout=subprocess.PIPE, shell=True, stderr=subprocess.STDOUT)
        out_buf = ret.stdout.read()
        self.logging.info(out_buf)

        self.models = BooleanModels("out/networks.csv")

        # TODO strategies
        #script = "caspo analyze --networks ./out/networks.csv  --midas {0} {1} --networks-stats".format(midasname.name, timepoint)

        #self.logging.info(script)

        #ret = subprocess.Popen("%s" % (script),
        #    stdout=subprocess.PIPE, shell=True, stderr=subprocess.STDOUT)
        #out_buf = ret.stdout.read()

        # ignore time zero in the computation so to agree with CNO, we divide
        # MSE by 2
        #self.logging.info(out_buf)
        #try:
        #    mse = float([x for x in out_buf.split("\n") if 'MSE' in x][0].split()[2])
        #except:
        #    mse = []

        # TODO simplify
        import pandas as pd
        self.models.scores = pd.read_csv("out/networks-mse-len.csv")['MSE']
        self.models.midas = self.midas

        #all_sizes = BooleanModels("out/networks-mse-len.csv").df['SIZE']

        self.logging.info("Best MSE: {0}".format(self.models.scores.min()/2.))

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
        from cno.core.results import BooleanResults
        br = BooleanResults()
        br.models = self.models
        br.hist_scores()

    def _get_scores_vs_model_size_df(self):
        df = pd.DataFrame(
            {'score': self.models.scores.values,
            'size': self.models.sizes.values})
        return df

    def scatter_scores_vs_model_size(self):
        df = self._get_scores_vs_model_size_df()
        viz.scatter_hist(df)

    def hist2d_scores_vs_model_size(self, bins=[10,10], cmap='gist_heat_r',fontsize=16):
        from cno.core.results import BooleanResults
        br = BooleanResults()
        br.models = self.models
        br.hist2d_scores_vs_model_size(cmap=cmap, bins=bins)

    def plot_model_number_vs_tolerance(self, normed=True, **kargs):
        """Plots number of models found as a function of the tolerance

        Tolerance is :math:`MSE_0 + fit`

        .. plot::
            :include-source:

            >>> from cno.boolean.asp import CASPO
            >>> from cno.datasets import cnodata
            >>> caspo = CASPO(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
            >>> caspo.optimise(fit=1, size=10)
            >>> a.plot_number_models_vs_tolerance()

        """
        from numpy import cumsum
        counter = Counter(self.models.scores)
        m = self.models.scores.min()
        counts = [counter[x] for x in sorted(counter.keys())]
        tolerances = [(x-m)/m*100 for x in sorted(counter.keys())]

        pylab.clf()
        marker = kargs.get('marker', 'o')
        linestyle = kargs.get('linestyle', '-')
        kargs['marker'] = marker
        kargs['linestyle'] = linestyle

        counts = cumsum(counts)
        if normed:
            counts /= max(counts)
        pylab.plot(tolerances, counts, **kargs)
        pylab.grid(True)
        pylab.xlabel("Tolerance in %")
        pylab.ylabel("Number of models")

    def __str__(self):
        """Prints information about the optimisation"""
        mses = self.models.scores.values
        #models = self.models[:]

        N = len(set(mses))
        txt = "There are %s different MSE found amongst %s models\n" % (N,len(mses))
        txt += "The minimum MSE is %s " % self.models.scores.min()
        return txt
        #Ngtts = len(self.family.gtts)
        #print("There are %s GTTS " % Ngtts)

    def summary(self):
        print(self)


    def get_exclusive_conjunctions(self):
        """Returns the exclusive conjunctions


        Some pairs of modules in models are mutually exclusive. You can
        then check whether replacing a module of each pair by the other
        has an effect or not on the MSE.


        """
        raise NotImplementedError
        conjunctions = list(self.family.combinatorics('exclusive'))
        for conjunction in conjunctions:

            ca = conjunction['conjunction_A'].__str__()
            cb = conjunction['conjunction_B'].__str__()
            fa = conjunction['frequency_A']
            fb = conjunction['frequency_B']
            print("CA: %s (frequency=%s) ; CB: %s (frequency=%s" % (ca, fa, cb, fb))

        return conjunctions

    def plot_gtts(self):
        """Plot sub-optimal models grouped by MSEs and GTTS

        Running the analysis with fit=0.1 like in the paper takes 8 hours...so
        don't be too greedy.

        :return: dictionary where keys are MSEs and values is a list of GTTS counting


        .. plot::
            :include-source:

            >>> from cellnopt.optimiser import ASPBool
            >>> from cellnopt.data import cnodata
            >>> a = ASPBool(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
            >>> a.run(fit=0.7, gtts=True)
            >>> res = a.plot_gtts()
            >>> res
            defaultdict(<type 'list'>, {0.0398: [2], 0.0486: [3, 2], 0.0426: [2], 0.0296: [1]})

        Here are the results for the ExtLiver data set and fit=0.1

        ::

            0.0499 [16]
            0.0507 [144, 16, 16]
            0.051 [518, 92, 92, 8, 8, 4, 2, 2]
            0.0519 [3126, 702, 702, 108, 108, 44, 34, 34, 16, 8, 8, 4, 4, 2, 2, 2, 2, 2, 2]
            0.0522 [92, 8, 8, 4]
            0.0523 [92, 8, 8]
            0.053 [702, 108, 108, 64, 8, 8, 8, 4, 4, 4, 4, 2, 2]
            0.0531 [112, 92, 16, 8, 8]
            0.0534 [52, 8, 4, 4]
            0.0539 [712, 144, 16, 8, 8, 4, 4]
            0.0542 [2090, 518, 426, 92, 60, 60, 36, 30, 30, 16, 8, 8, 8, 4, 4, 4, 2, 2, 2, 2, 2, 2]
            0.0543 [16]
            0.0546 [4]


        """
        raise NotImplementedError
        # A given GTT has a unique MSE. Several GTT may have the same MSE
        # We want to plot the figure 5 of ASP/CASO paper so we want to extract
        # for a given MSE, all the GTTs and figure out how many models are found
        # in each of them.

        # get all GTTs and get their respective MSEs.
        mses, Ngtts = self._get_mse_Ngtts()

        # what are the unique MSEs ?
        # round values to agree with figure 4
        unique_mses = set([round(x,4) for x in mses])
        group = defaultdict(list)

        # Now, for each unique MSE, retrieve the different GTT and their number
        for this_mse in unique_mses:
            for i, mse in enumerate(mses):
                if round(mse, 4) == this_mse:
                    print( this_mse, Ngtts[i])
                    group[this_mse].append(Ngtts[i])
            print("----")

        for this in sorted(unique_mses):
            print(this, sorted(group[this], reverse=True))

        pylab.clf()
        N = 0
        icolor = 0
        colors = ['yellow','k','b','r','g', 'orange', '', '', '', '']

        M = max(Ngtts) * 1.1

        for m in sorted(unique_mses):
            print(m, sorted(group[m], reverse=True))
            pylab.bar(range(N, N+len(group[m])), sorted(group[m], reverse=True,),
                   width=0.8 , color=colors[icolor%6], label="%s" % m) # %6 = number of colors
            pylab.plot([N,N],[0, M], 'k--')

            icolor += 1
            print(N, N+len(group[m]))
            N += len(group[m])
        pylab.grid(axis="y")
        pylab.legend(loc="upper left", ncol=2)
        pylab.xticks([],[])
        pylab.ylabel("#")
        pylab.title("""Distribution of sub-optimal models
        Models are ordered in groups (colored) from left to right
        according to their MSEs. Each bar within a group represents a GTT""")
        return group

    def _get_mse_Ngtts(self):
        if self.gtts:
            Ngtts = [ len(x) for x in list(self.family.gtts)]
            mses = [x.mse(self.dataset)/self.normalisation for x in list(self.family.gtts)]
            return (mses, Ngtts)
    def _get_mse_gtt(self):
        if self.gtts:
            Ngtts = [ len(x) for x in list(self.family.gtts)]














