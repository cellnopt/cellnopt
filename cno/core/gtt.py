# -*- python -*-
#
#  This file is part of the CNO package
#
#  Copyright (c) 2012-2013 - EMBL-EBI
#
#  File author(s): Thomas Cokelaer (cokelaer@ebi.ac.uk)
#
#  Distributed under the GLPv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: github.com/cellnopt/cellnopt
#
##############################################################################
import numpy as np
import pylab
import pandas as pd
from cno.core import CNORBase
from easydev import TempFile


__all__ = ['GTTBool']


class GTTBool(object):
    """API could be simplified.

    ::


        # here there are duplicate, should be called automatically ?
        s = steady.Steady(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
        s.preprocessing()
        s.optimise()
        s.results.models.drop_duplicates()


        g = gtt.GTTBool(s)
        g.analyse()
        # here indices are same s in models, which may not be contiguous
        g.truth_tables[189]


    """

    def __init__(self, simulator):
        self.simulator = simulator # do not touch

    def analyse(self):
        models = self.simulator.results.models
        self.truth_tables = {}
        from easydev import progress_bar
        pb = progress_bar(len(models.df))
        for i, index in enumerate(models.df.index):
            reactions = models.df.ix[index][models.df.ix[index]==1]
            reactions = list(reactions.index)
            self.simulator.simulate(reactions=reactions)
            tt = self.simulator.simulated[self.simulator.time].flatten()
            self.truth_tables[index] = tt
            pb.animate(i+1)







class GTTBoolOld(CNORBase):
    """Works for boolean steady state (1 time point)

    ::

        from cno import *
        c = cnorbool.CNORbool(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"), verboseR=False)
        c.optimise(reltol=0.5)
        c.optimise(reltol=0.5)
        g   = gtt.GTTBool(c._model, c.data, c.models, c.results.scores)
        d = g.get_gtt()
    """
    def __init__(self, model, data, models, scores, verboseR=False):
        """
        Note that once, grouped, the scores should be identical albeit the
        model size

        [scores[i] for i in grouped.groups.values()[10]]


        :param model: a instance of :class:`CNOGraph`
        :param data: an instance of :class:`XMIDAS`
        :param models: an instance of compatible :class:`Models`
        :param scores: the scores of each model.

        """
        CNORBase.__init__(self, verboseR)
        self.models = models
        self.scores = scores
        self.model = model
        self.data = data # a MIDAS file

    def _init(self):
        """

        preprocessing in R and Python may lead to different labels (not logic)
        So, we need to reprocess unfortunately.

        """
        fhmodel = TempFile()
        fhdata = TempFile()
        self.model.to_sif(fhmodel.name)
        self.data.to_midas(fhdata.name)
        self.session.run("library(CellNOptR)")
        self.session.run('model=readSIF("%s")' % fhmodel.name)
        self.session.run('cnolist=CNOlist("%s")' % fhdata.name)

    def _get_sim(self, bs):
        self.session.bs1 = bs
        script = """
        png()
        output = cutAndPlot(cnolist, model, list(bs1), plotPDF=F)
        dev.off()
        """
        self.session.run(script)
        res = self.session.output['simResults'][0]
        res = list(res['t0'].flatten() ) + list(res['t1'].flatten())
        return res

    def compute_gtts(self):
        print("init R library")
        self._init()
        N = len(self.models)
        from easydev import progress_bar
        b = progress_bar(N)
        d = {}
        for i in range(0, N):
            res = np.array(self._get_sim(self.models.df.ix[i].values))
            b.animate(i)
            d[i] = res

        df = pd.DataFrame(d).transpose()
        grouped = df.groupby(list(df.columns))
        pylab.hist([len(this) for this in grouped.groups.values()], 100)
        res =  {'df':df, 'simulation': d, 'grouped':grouped}#
        self.gtts = res
        return self.gtts

    def plot_gtt(self, bins=100):
        grouped = self.gtts
        pylab.hist([len(this) for this in grouped.groups.values()], bins)

    def __len__(self):
        return len(self.scores)


class GTTDendrogram(object):

    def __init__(self, asp):

        self._asp = asp # reference that must not be touch 
        self.setup = Setup(
            self._asp.cnograph.stimuli, 
            self._asp.cnograph.inhibitors, 
            self._asp.cnograph.signals)
    
        # Set the experiments 
        self.experiments = []
        ac = self._asp.cnograph
        for d in list(self._asp.dataset):

            inputs = [(x[0],x[1]) if x[0] not in ac.inhibitors else (x[0]+"i", x[1]) 
                    for x in d.inputs.iteritems()]
            self.experiments.append(dict(inputs))

        # Set the outputs to be looked at
        self.outputs = sorted(self._asp.dataset[0].outputs[0].keys())

    def _get_asp(self):
        return self._asp
    asp = property(_get_asp)

    def _get_readouts_for_one_gtt(self, gtt):
        """# create experiments given setup for a given GTT
        # keep only experiments that we are interested in (ie present in dataset
        # experiments."""
        tokeep = []
        names = self.experiments[0].keys()

        for i,r in enumerate(gtt.to_matrix(self.setup)):
            thisexp = dict([x for x in r.iteritems() if x[0] in names])
            #if thisexp in experiments:

            if thisexp in self.experiments:
                #print i,
                #res = [thisexp[x] for x in names] + [r[x] for x in outputs]
                # we need only the outputs since inputs/experiments are always
                # the same
                res = [r[x] for x in self.outputs]
                tokeep.append(res)

        m = matrix(tokeep).flatten()
        return m


    # get matrix for the first gtt
    def get_readouts(self):
        print("Getting all readouts for all %s GTTS" % len(self._asp.family.gtts))
        print("Using the following outputs %s" % self.outputs)

        allres = []
        for i, gtt in enumerate(list(self._asp.family.gtts)):
            print "looking at GTT %s " % i
            m = self._get_readouts_for_one_gtt(gtt)
            allres.append(m)
        N = len(self._asp.family.gtts)
        M = len(self.outputs) * len(self.experiments)


        m = array(allres).reshape(N, M)

        return m


    def get_distance_matrix(self, *args, **kargs):
        allres = self.get_readouts()

        # now compute the dist between the vectors
        dm = pdist(allres, *args, **kargs)
        return dm

    def get_linkage(self, dm):
        linkageMatrix = linkage(dm)
        return linkageMatrix

    def dendrogram(self, linkage=None, metric="euclidean", **kargs):
        """

            res = d.dendrogram(metric="hamming", 
                        leaf_label_func=lambda x:d.get_mses()[x], 
                        leaf_rotation=90)

                        leaf_label_func=lambda x:d.get_GTTS_length[x], 
        """
        if linkage == None:
            dm = self.get_distance_matrix(metric=metric)
            linkage = self.get_linkage(dm)

        from pylab import clf
        clf()
        res = dendrogram(linkage, show_leaf_counts=True, leaf_font_size=16, **kargs)
        return res

    def get_mses(self):
        mses = [round(x.mse(self._asp.dataset), 4) for x in self._asp.family.gtts]
        return mses
    def get_GTTs_length(self):
        N = [len(gtt.models) for gtt in self._asp.family.gtts]
        return N

    
