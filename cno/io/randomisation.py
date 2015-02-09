
from cno.io import xcnograph


class RandomGraph(object):
    """

    c = RandomGraph()
    c.run(10, 200)
    c.parameters['nExtraNode'] = 20

    """
    def __init__(self):
        self.commons = []
        self.parameters = {'nStimuli':4, 'nSignals':4, 'extraNode':10, 'nTrans':4}    

    def create_graph(self, fraction=0.7):
        self.c = xcnograph.XCNOGraph()
        self.c.random_cnograph(nStimuli=self.parameters['nStimuli'], 
                nSignals=self.parameters['nSignals'], 
                fraction_activation=fraction, nTranscript=self.parameters['nTrans'], 
                nExtraNode=self.parameters['extraNode'])
        self.Nneg = [x[2]['link'] for x in self.c.edges(data=True)].count('-')
        self.indeg = self.c.in_degree().copy()
        self.outdeg = self.c.out_degree().copy()
        self.deg = self.c.degree().copy()
        self.ori = self.c.copy()
        self.c = self.c.copy()
        self.Nedges = len(self.ori.edges())

    def swap(self, N=100):
        common = []
        for i in range(0, N):
            self.c.swap_edges(1, self_loop=False)
            #c.plot(rank_method='cno')

            assert self.Nneg == [x[2]['link'] for x in self.c.edges(data=True)].count('-')
            assert self.deg == self.c.degree()
            assert self.indeg == self.c.in_degree()
            assert self.outdeg == self.c.out_degree()
            n = len([e for e in self.c.edges() if e in self.ori.edges()])/float(self.Nedges)
            common.append(n)
        self.commons.append(common)

    def plot(self):
        from pylab import plot, show, grid
        for common in self.commons:
            plot(common)    

        import pandas as pd
        df = pd.DataFrame(self.commons)
        df = plot(df.mean(), lw=2, color='r')

    def run(self, N=10, Nswaps=100, verbose=True):
        self.commons = []
        for i in range(0,10):
            if verbose:
                print(i)
            self.create_graph()
            self.swap(N=Nswaps)
        self.plot()
        import pylab
        pylab.grid()




