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
from easydev import precision

and_symbol = "^"


__all = ["Models", "BooleanModels", "ContinousModels", 
    "DTModels", "FuzzyModels", "CompareModels"]


class Models(object):
    """Data structure to store models.

    Models are stored in dataframes. Columns will hold the reactions.


    """
    def __init__(self, data, reacID=None, index_col=None, verbose=True):
        """.. rubric:: constructor

        :param data: could be a string to read a file (CSV). The CSV header 
            should contain the reaction name. First column is not expected 
            to be found. Note, however that :param:`index_col` may be used.
            The input can also be a dataframe (column names being the reactions)
            set to 0/1. The input can also be an instance of :class:`Models`.
        :param list reacID: if provided, columns are renamed using this list
        :param index_col:
        :param bool verbose:

        Reaction names may contain a symbol indicating the logical ANDs. This
        should be "^" character.

        """
        self.verbose = verbose
        # FIXME interpret the first columns automatically ?
        if isinstance(data, str):
            self.filename = data
            self.df = pd.read_csv(self.filename, index_col=index_col)
            # FIXIME What is this
            if reacID:
                reacID = pd.read_csv(reacID)
                self.df.columns = reacID.ix[:,0]

            if 'Score' in self.df.columns:
                self.scores = self.df.Score
                del self.df['Score']

            if 'score' in self.df.columns:
                self.scores = self.df.score
                del self.df['score']
        elif isinstance(data, pd.DataFrame):
            self.df = data.copy()
        elif isinstance(data, Models):
            self.df = data.df.copy()
        else:
            from cno import CNOError
            raise CNOError("input data not understood. Could be a filename, a dataframe or a Models instance")

        if hasattr(data, 'scores'):
            self.scores = getattr(data, 'scores')

        # TODO: In a reaction from cnograph, they should be not ORs, just simple
        # reactions and ANDS (e.g., A^B=C). If "A+B=C" is found, this is coming
        # from CellNOptR, ,which has a different conventions. So, we replace
        # all + by "^" !! Do we want a warning ?
        for reaction in self.df.columns:
            count = 0
            if "+" in reaction:
                # todo: use logging

                if self.verbose and count == 0:
                    print("Warning in Models. found a + sign... in %s. Interepreted as ^" % reaction)
                    count = 1
            def convert(x):
                from cno import Reaction
                r = Reaction(x)
                r.sort()
                name = r.name
                name = name.replace("+", "^")
                return name
            self.df.columns = [convert(x) for x in self.df.columns]

        # we also reorder alphabetically the species in the and reactions

        # keep this import here to avoid cycling imports
        from cno.io.cnograph import CNOGraph
        from cno.io import Reaction
        self.cnograph = CNOGraph()
        non_reactions = []
        for this in self.df.columns:
            try:
                reac = Reaction(str(this))
                self.cnograph.add_reaction(str(this))
            except:
                if self.verbose:
                    print('Skipping column %s (not valid reaction ?)' % this)
                non_reactions.append(this)
        #self.non_reactions = non_reactions
        #self.df_non_reactions = self.df[non_reactions].copy()

    def drop_scores_above(self, tolerance=None):
        max_score = self.scores.min() * (1+tolerance)
        index = self.df.ix[self.scores<=max_score].index

        self.df = self.df.ix[index]
        self.scores = self.scores.ix[index]

    def get_average_model(self, max_score=None):

        """Returns the average model (on each reaction)"""
        if max_score is None:
            return self.df.mean(axis=0)
        else:
            #filter scores below some vlues
            N = float(sum(self.scores<=max_score))
            print('Keeping %s percent of the models' % str( N /len(self.scores)*100.))

            return self.df.ix[self.scores<=max_score].mean(axis=0)

    def to_csv(self, filename, index=False):
        """Exports the dataframe to a CSV file"""
        self.df['score'] = self.scores.values
        self.df.to_csv(filename, index=False)
        del self.df['score']

    def to_sif(self, filename=None):
        """Exports 2 SIF using the "and" convention

        can read the results with CellNOptR for instance::

            library(CellNOptR)
            plotModel(readSIF("test.sif"))

        """
        return self.cnograph.to_sif(filename)

    def __eq__(self, other):
        if len(self.df) != len(other.df):
            return False
        df1 = self.df.copy()
        df2 = other.df.copy()
        if all(df1.columns != df2.columns):
            return False
        # make sure the columns are ordered similarly
        df2 = df2[df1.columns]
        return all(df1.sort() == df2.sort())

    def __len__(self):
        return len(self.df)


class FuzzyModels(Models):
    def __init__(self, data, reacID=None, index_col=None):
        super(FuzzyModels, self).__init__(data, reacID, index_col)

    def copy(self):
        return FuzzyModels(self)


class BooleanModels(Models):
    """Class to read and plot models as exported by CASPO or CellNOptR


    Models contains dataframe with reactions as columns and models as rows.
    For each reaction, we can then obtain the average paramters for a reaction.
    In a boolean case, a Model stores a value made of 0/1

    scores may be available. No sizes are stored. Sizes could be extracted easily
    as sum over rows.

    ::

        >>> from cno.core.models import Models
        >>> m = Models()
        >>> m.plot() # average model, whcih can be obtained with  m.get_average_model()
        >>> m.plot(model_number=0)  # indices are m.df.index
        >>> m.plot(model_number=0)  # indices are m.df.index

    .. note:: One difficulty is the way ANDs are coded in different software. In CASPO,
        the AND gate is coded as "A+B=C". Note that internally we use ^ especially
        in CNOGraph. Then, an AND edge is splitted in sub edges. so, A+B=C is made
        of 3 edges A -> A+B=C , B -> A+B=C and A+B=C -> C. This explains the wierd
        code in :meth:`cno.io.cnograph.plot`.


    - plots average models with edges on/off
    - plot of errobars on edges sorted by average presence
    - plots heatmap of the models
    """
    def __init__(self, data, reacID=None, index_col=None):
        """
        if you have a first column, whihc is not a reaction, set index_col to 0

        .. todo:: values are 0/1 since we have bit strings but could be anything in other
            formalisms (e.g., ODE) how to handle those cases ?

        :param dta: a filename with columns as the reacitons and rowss as
            parameters for each reactions. Each row is therefore a model.
        """
        super(BooleanModels, self).__init__(data, reacID, index_col)

    def get_cv_model(self):
        """Returns the average coefficient of variation on each reaction"""
        res = self.df.std(axis=0)/self.df.mean(axis=0)
        res = res.fillna(0)
        return res

    def compute_average(self, model_number=None, tolerance=None):
        """Compute the average and update the cnograph accordingly

        :param int model_number: model_number as shown by :attr:`df.index`
            if not provided, the average is taken
        """
        if model_number is None and tolerance is None:
            model = self.get_average_model()
        elif model_number == 'cv':
            model = self.get_cv_model()
        elif tolerance is not None:
            model = self.get_average_model(max_score = self.scores.min() * (1.+tolerance))
            if len(model) == 0:
                raise ValueError('No model found within that tolerance')
        else:
            model = self.df.ix[model_number]

        # This is to set the average and label and penwidth
        # TODO: could be simplified using Reaction ?
        for edge in self.cnograph.edges(data=True):
            link = edge[2]['link']
            if and_symbol not in edge[0] and and_symbol not in edge[1]:
                if link == "-" :
                    name = "!" + edge[0] + "=" + edge[1]
                else:
                    name = edge[0] + "=" + edge[1]
                value = model[name]
            elif and_symbol in edge[0]:
                value = model[edge[0]]
            elif and_symbol in edge[1]:
                value = model[edge[1]]
            else:
                raise ValueError()
            self.cnograph.edge[edge[0]][edge[1]]["label"] = precision(value)
            self.cnograph.edge[edge[0]][edge[1]]["average"] = precision(value)

            # if values are between 0 and 1
            M = float(model.max())
            self.cnograph.edge[edge[0]][edge[1]]["penwidth"] = precision(value, 2) * 5/M

    def plot(self, model_number=None, cmap='gist_heat_r',
            colorbar=True, tolerance=None, filename=None, **kargs):
        """Plot the average model"""
        self.compute_average(model_number=model_number, tolerance=tolerance)
        self.cnograph.plot(edge_attribute="average", cmap=cmap,
                colorbar=colorbar, filename=filename, **kargs)

    def errorbar(self, tolerance=1e8, errorbar=True):
        """Plot the average presence of reactions over all models"""

        try:
            df = self.df.ix[self.scores<=self.scores.min()*(1+tolerance)]
        except:
            df = self.df[(self.scores<=self.scores.min()*(1+tolerance)).values]

        mu = df.mean()
        mu.sort(inplace=True)
        sigma = df.std()
        pylab.clf()
        X = range(0,len(mu.index))
        if errorbar is True:
            errorbar = 1
        else:
            errorbar = 0
        pylab.errorbar(X, mu.values, yerr=sigma.ix[mu.index].values*errorbar,
                       marker='x', color='r', lw=0, elinewidth=2, ecolor='b')
        pylab.xticks(X, mu.index, rotation=90)
        pylab.title('')
        pylab.grid()
        pylab.ylim([-0.1, 1.1])
        #pylab.xlim([-0.5, len(X)+.5])
        pylab.tight_layout()
        return df

    def heatmap(self, num=1, transpose=False, cmap='gist_heat_r', heatmap_attr={}):
        """    """
        #df = self.get_average_models()
        from biokit.viz.heatmap import Heatmap
        if transpose:
            df = self.df.transpose()
        else:
            df = self.df
        h = Heatmap(df)
        h.plot(cmap=cmap,num=num, **heatmap_attr)
        return h

    def __add__(self, other):
        import pandas as pd
        df = pd.concat([self.df, other.df])
        df.drop_duplicates(inplace=True)
        return Models(df)

    def __str__(self):
        txt = "Models contains {0} rows".format(len(self))
        return txt

    def copy(self):
        return BooleanModels(self)

    def _get_sizes(self):
        return self.df.sum(axis=1)    
    sizes = property(_get_sizes)

    def drop_duplicates(self):
        self.df['score'] = self.scores
        self.df.drop_duplicates(inplace=True)
        self.scores = self.df['score']
        del self.df['score']

    def get_main_reactions(self, threshold=0.5):
        reactions = list(self.df.columns[self.df.mean() > threshold])
        reactions = [x.replace('+','^') for x in reactions]
        return reactions

    def get_consensus_model(self, threshold=0.5):
        df = self.df.ix[self.scores<=self.scores.min()*(1.)]
        reactions = list(df.mean()[df.mean() > threshold].index)
        return reactions

    def get_jaccard(self, progress=True):
        import sklearn.metrics
        N = len(self.df)
        J = np.zeros((N,N))
        from easydev import progress_bar
        pb = progress_bar(N)
        for ic, i in enumerate(self.df.index):
            for jc, j in enumerate(self.df.index):
                J[ic][jc] = sklearn.metrics.jaccard_similarity_score(self.df.ix[i], self.df.ix[j])
            pb.animate(1+ic)
        return J





class DTModels(BooleanModels):
    def __init__(self, data, reacID=None, index_col=None):
        super(DTModels, self).__init__(data, reacID, index_col)

    def copy(self):
        return DTModels(self)


class CompareTwoModels(object):
    """

    """
    def __init__(self, m1, m2):
        """

        :param m1: first model as a Pandas time series e.g. row of BooleanModels
        :param m2: first model as a Pandas time series e.g. row of BooleanModels
        :return:

        from a models, m1 = pd.TimeSeries(models.df.ix[0], dtype=int)
        m2 = pd.TimeSeries(models.df.ix[1], dtype=int)
        """
        
        self.m1 = m1
        self.m2 = m2

        assert all(self.m1.index == self.m2.index) == True
        self.midas = None

    def get_intersection(self):
        return self.m1[np.logical_and(self.m1, self.m2)]

    def get_union(self):
        return self.m1[np.logical_or(self.m1 , self.m2)]

    def get_both(self):
        return self.get_intersection()

    def get_m1_only(self):
        return self.m1[np.logical_and(self.m1==1, self.m2==0)]

    def get_m2_only(self):
         return self.m2[np.logical_and(self.m1==0, self.m2==1)]

    def get_both_off(self):
         return self.m2[np.logical_and(self.m1==0, self.m2==0)]

    def plot_multigraph(self, cmap='jet'):
        from cno.io.multigraph import  CNOGraphMultiEdges
        #from cno import CNOGraph
        from cno import Reaction

        c = CNOGraphMultiEdges()
        c.midas = self.midas

        for reaction in self.get_both().index:
            r = Reaction(reaction)
            r.sort()
            for edge, link in c.reac2edges(r.name):
                c.add_edge(edge[0], edge[1], link=link, edgecolor=.1, color='black', penwidth=6, label='both')

        for reaction in self.get_m1_only().index:
            r = Reaction(reaction)
            r.sort()
            for edge, link in c.reac2edges(r.name):
                c.add_edge(edge[0], edge[1], link=link, edgecolor=.3, label='m1', color='red', penwidth=3)

        for reaction in self.get_m2_only().index:
            r = Reaction(reaction)
            r.sort()
            for edge, link in c.reac2edges(r.name):
                c.add_edge(edge[0], edge[1], link=link, edgecolor=.5, label='m2', color='green', penwidth=3)

        for reaction in self.get_both_off().index:
            r = Reaction(reaction)
            r.sort()
            for edge, link in c.reac2edges(r.name):
                c.add_edge(edge[0], edge[1], link=link, edgecolor=.9, label='', arrowsize=0, color='gray', penwidth=0)


        #c.plot(edge_attribute='edgecolor', cmap=cmap)
        c.plot()
        return c




class MultiModels(object):
    pass




