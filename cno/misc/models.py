import pylab
import pandas as pd
from easydev import precision

and_symbol = "^"



class Models(object):
    """Class to read and plot models as exported by CASPO or CellNOptR


    Models contains dataframe with reactions as columns and models as rows.
    For each reaction, we can then obtain the average paramters for a reaction.
    In a boolean case, a Model stores a value made of 0/1

    No scores are stored. No sizes are stored. Sizes could be extracted easily 
    as sum over rows.

    ::

        >>> import models
        >>> m = models.Models()
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
        # FIXME interpret the first columns automatically ?
        if isinstance(data, str):
            self.filename = data
            self.df = pd.read_csv(self.filename, index_col=index_col)
            if reacID:
                reacID = pd.read_csv(reacID)
                self.df.columns = reacID.ix[:,0]
        elif isinstance(data, pd.DataFrame):
            self.df = data.copy()
        elif isinstance(data, Models):
            self.df = data.df.copy()
        else:
            raise CNOError("input data not understood. Could be a filename, a dataframe or a Models instance")

        # TODO: In a reaction from cnograph, they should be not ORs, just simple
        # reactions and ANDS (e.g., A^B=C). If "A+B=C" is found, this is coming
        # from CellNOptR, ,which has a different conventions. So, we replace
        # all + by "^" !! Do we want a warning ?
        self.df.columns = [x.replace("+", "^") for x in self.df.columns]

        # keep this import here to avoid cycling imports
        from cno.io.cnograph import CNOGraph
        self.cnograph = CNOGraph()
        for this in self.df.columns:
            self.cnograph.add_reaction(this)

    def get_average_model(self):
        """Returns the average model (on each reaction)"""
        return self.df.mean(axis=0)

    def get_cv_model(self):
        """Returns the average coefficient of variation on each reaction"""
        res = self.df.std(axis=0)/self.df.mean(axis=0)
        res = res.fillna(0)
        return res

    def compute_average(self, model_number=None):
        """Compute the average and update the cnograph accordingly

        :param int model_number: model_number as shown by :attr:`df.index`
            if not provided, the average is taken
        """
        if model_number is None:
            model = self.get_average_model()
        elif model_number == 'cv':
            model = self.get_cv_model()
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
            colorbar=True, *args, **kargs):
        """Plot the average model"""
        self.compute_average(model_number=model_number)
        self.cnograph.plot(edge_attribute="average", cmap=cmap, 
                colorbar=colorbar,**kargs)

    def to_csv(self, filename):
        self.df.to_csv(filename)

    def to_sif(self, filename=None):
        """Exports 2 SIF using the "and" convention

        can read the results with CellNOptR for instance

            >>> library(CellNOptR)
            >>> plotModel(readSIF("test.sif"))
        """
        return self.cnograph.to_sif(filename)

    def errorbar(self):
        """Plot the average presence of reactions over all models"""
        mu = self.df.mean()
        mu.sort(inplace=True)
        sigma = self.df.std()
        pylab.clf()
        X = range(0,len(mu.index))
        pylab.errorbar(X, mu.values, yerr=sigma.ix[mu.index].values,
                       marker='x', color='r', lw=0, elinewidth=2, ecolor='b')
        pylab.xticks(X, mu.index, rotation=90)
        pylab.title('')
        pylab.grid()
        pylab.ylim([-0.1, 1.1])
        pylab.xlim([-0.5, len(X)+.5])
        pylab.tight_layout()

    def heatmap(self, num=1, transpose=False, cmap='gist_heat_r', heatmap_attr={}):
        """

        .. plot::
            :include-source:

            from corda import *
            m = Models(fit=0, factor=1)
            m.heatmap()
        """
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

    def __str__(self):
        txt = "Models contains {0} rows".format(len(self))
        return txt

    def copy(self):
        return Models(self)

