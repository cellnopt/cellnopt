import pylab
import pandas as pd
from easydev import precision


class Models(object):
    """Class to read and plot models as exported by CASPO

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


    """
    def __init__(self, data, reacID=None, index_col=None):
        """
        if you have a first column, whihc is not a reaction, set index_col to 0

        .. todo:: values are 0/1 since we have bit strings but could be anything in other
        formalisms (e.g., ODE) how to handle those cases ?
        """
        # FIXME interpret the first columns
        if isinstance(data, str):
            self.filename = data
            self.df = pd.read_csv(self.filename, index_col=index_col)
            if reacID:
                reacID = pd.read_csv(reacID)
                self.df.columns = reacID.ix[:,0]
        elif isinstance(data, pd.DataFrame):
            self.df = data.copy()

        self.df.columns = [x.replace("+", "^") for x in self.df.columns]

        # keep this import here to avoid cycling imports
        from cno.io.cnograph import CNOGraph
        self.cnograph = CNOGraph()
        try:
            for this in self.df.columns:
                self.cnograph.add_reaction(this)
        except:
            print("could not interpret some reactions. cnograph may not be valid")
            pass

    def get_average_model(self):
        """Returns the average model (on each reaction)"""
        return self.df.mean(axis=0)

    def get_cv_model(self):
        """Returns the average coefficient of variation on each reaction"""
        res = self.df.std(axis=0)/self.df.mean(axis=0)
        res = res.fillna(0)
        return res

    def compute_average(self, model_number=None, *args, **kargs):
        """Compute the average and update the cnograph accordingly

        :param int model_number: model_number as shown by :attr:`df.index`
            if not provided, the average is taken
        """
        if model_number==None:
            model = self.get_average_model()
        elif model_number == 'cv':
            model = self.get_cv_model()
        else:
            model = self.df.ix[model_number]

        for edge in self.cnograph.edges(data=True):
            link = edge[2]['link']
            if "^" not in edge[0] and "^" not in edge[1]:
                if link=="-":
                    name = "!"+edge[0]+"="+edge[1]
                else:
                    name = edge[0]+"="+edge[1]
                value = model[name]
            elif "^" in edge[0]:
                value = model[edge[0]]
            elif "^" in edge[1]:
                value = model[edge[1]]
            else:
                raise ValueError()
            self.cnograph.edge[edge[0]][edge[1]]["label"] = precision(value)
            self.cnograph.edge[edge[0]][edge[1]]["average"] = precision(value)

            # if values are between 0 and 1
            M = float(model.max())
            self.cnograph.edge[edge[0]][edge[1]]["penwidth"] = precision(value,2)*100/M/10.

    def plot(self, model_number=None, *args, **kargs):
        """Plot the average model

        """
        self.compute_average(model_number=model_number, *args, **kargs)
        self.cnograph.plot(edge_attribute="average", **kargs)

    def to_sif(self, filename=None):
        """Exports 2 SIF using the "and" convention

        can read the results with CellNOptR for instance

            >>> library(CellNOptR)
            >>> plotModel(readSIF("test.sif"))
        """
        return self.cnograph.to_sif(filename)

    def errorbar(self):
        """Plot the average presence of reactions over all models


        """
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

