import pylab
import pandas as pd
from easydev import precision


class Models(object):
    """Class to read and plot models as exported by CASPO

        >>> import models
        >>> m = models.Models()
        >>> m.plotdot() # average model, whcih can be obtained with  m.get_average_model()
        >>> m.plotdot(model_number=0)  # indices are m.df.index
        >>> m.plotdot(model_number=0)  # indices are m.df.index

    .. note:: One difficulty is the way ANDs are coded in different software. In CASPO,
        the AND gate is coded as "A+B=C". Note that internally we use ^ especially
        in CNOGraph. Then, an AND edge is splitted in sub edges. so, A+B=C is made
        of 3 edges A -> A+B=C , B -> A+B=C and A+B=C -> C. This explains the wierd
        code in :meth:`plotdot`.

    ::

    #Manually:

    from cellnopt.core import *
    c = CNOGraph("TCellcombiMS2.sif")
    import pandas as pd
    pd.read_csv("reacIDs.csv")
    networks.columns = reacIDs.values
    reacIDs = pd.read_csv("reacIDs.csv")
    networks = pd.read_csv("networks.csv")
    networks.columns = reacIDs.ix[:,0]

    # or using cellnopt.core.models

    from cellnopt.core import models
    m = models.Models("networks.csv", "reacIDs.csv")
    # a cnograph is built inside m.cnograph, which should be identical to the input SIF file
    m.plotdot()
    m.edges()
    for e in m.cnograph.edges(data=True):
        e[2]['label'] = int(e[2]['label']*100)/100.
    for e in m.cnograph.edges(data=True):
        e[2]['penwidth'] = e[2]['label']
    for e in m.cnograph.edges(data=True):
        e[2]['average'] = int(e[2]['average']*100)
    for e in m.cnograph.edges(data=True):
        e[2]['label'] = e[2]['average']

    for e in m.cnograph.edges(data=True):
        e[2]['penwidth'] = log2(1+ float(e[2]['average']))
    for e in m.cnograph.edges(data=True):
        e[2]['average'] = e[2]['average']/100.

    for e in m.cnograph.edges(data=True):
        e[2]['penwidth'] = log2(1+ float(e[2]['label']))
    for e in m.cnograph.edges(data=True):
        if e[2]['average']<0.02:
            m.cnograph.remove_edge(e[0], e[1])
    for n in m.cnograph.nodes():
        if m.cnograph.degree(n)==0:
            m.cnograph.remove_node(n)

    m.cnograph.plot(edge_attribute="average", cmap="gray", edge_attribute_labels=True)
    m.cnograph.midas.remove_stimuli("BN201")
    m.cnograph.midas.remove_stimuli("gilenya")
    m.cnograph.midas.remove_stimuli("IFNG")
    m.cnograph.midas.remove_stimuli("Teriflunomide")
    m.cnograph.midas.remove_stimuli("INS")
    m.cnograph.midas.remove_stimuli("S1P1")
    m.cnograph.midas.remove_stimuli("vitD3")
    m.cnograph.midas.remove_stimuli("DMF")
    for e in m.cnograph.edges(data=True):
        e[2]['penwidth'] = log2(1+ float(e[2]['label']))
    m.cnograph.plot(edge_attribute="average", cmap="copper", edge_attribute_labels=True)

    # remove edges with weight < 0.05. We need to get rid of nodes that are not connected anymore
    for e in m.cnograph.edges(data=True):
        if e[2]['average']<0.05:
           m.cnograph.remove_edge(e[0], e[1])
    m.cnograph.plot(edge_attribute="average", cmap="copper", edge_attribute_labels=True)
    for n in m.cnograph.nodes():
        if m.cnograph.degree(n)==0:
            m.cnograph.remove_node(n)
    m.cnograph.midas.remove_stimuli("antiCD3")
    m.cnograph.midas.remove_stimuli("BDNF")
    m.cnograph.midas.remove_stimuli("H202")
    m.cnograph.midas.remove_stimuli("H2O2")
    for e in m.cnograph.edges(data=True):
        e[2]['penwidth'] = log2(1+ float(e[2]['label']))
    m.cnograph.plot(edge_attribute="average", cmap="copper", edge_attribute_labels=True)


    """
    def __init__(self, data, reacID=None, index_col=None):
        """
        if you have a first column, whihc is not a reaction, set index_col to 0
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
        #if any(self.df>1):
        #    self.df /= 100.
        #self.df = self.df.ix[:,0:450]
        # cellnoptR uses + for ANDs
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
        """Returns the average model"""
        return self.df.mean(axis=0)


    def get_cv_model(self):
        """Returns the average model"""
        res = self.df.std(axis=0)/self.df.mean(axis=0)
        res = res.fillna(0)
        return res


    def compute_average(self, model_number=None, *args, **kargs):
        """

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
        self.compute_average(model_number=model_number, *args, **kargs)
        self.cnograph.plotdot(edge_attribute="average", **kargs)

    def export2sif(self, filename):
        """Exports 2 SIF using the "and" convention

        can read the results with CellNOptR for instance

            >>> library(CellNOptR)
            >>> plotModel(readSIF("test.sif"))
        """
        self.cnograph.export2sif(filename)

    def errorbar(self):
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

