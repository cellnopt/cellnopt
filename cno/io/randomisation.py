
from cno import *
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

        import pylab
        import pandas as pd
        N = len(self.commons[0])
        df = pd.DataFrame(self.commons)
        #df = plot(df.mean(), lw=2, color='r')

        pylab.clf()
        if N>200:
            pylab.errorbar(range(0,N)[::int(N/20)], df.mean()[::int(N/20)], df.std()[::int(N/20)])
        else:
            pylab.errorbar(range(0,N), df.mean(), df.std())
        pylab.plot(df.mean(), lw=2, color='r')

    def run(self, N=10, Nswaps=100, verbose=True):
        self.commons = []
        for i in range(0, N):
            if verbose:
                print("Creating graph %s " % i)
            self.create_graph()
            if verbose:
                print('Now the swapping')
            self.swap(N=Nswaps)
        self.plot()
        import pylab
        pylab.grid()


class Swaper(object):
    """Function to test the randomisation of networks.


    s = Swaper()
    s.plot_average_distance(10,100)

    """
    def __init__(self):

        self.mode = 'toy'
        self.init()

    def init(self):
        if self.mode == 'toy':
            self.graph = CNOGraph(cnodata("PKN-ToyMMB.sif"), 
                cnodata("MD-ToyMMB.csv"))
        else:
            self.graph = CNOGraph(cnodata("PKN-ExtLiverPCB.sif"), 
                cnodata("MD-ExtLiverPCB.csv"))
        self.edges = self.graph.edges()

    def swap(self, nswap=1, inplace=False):
        if inplace is True:
            self.graph.swap_edges(nswap, inplace=inplace)
        else:
            G, status = self.graph.swap_edges(nswap, inplace=inplace)
            return G

    def get_distance(self, G):

        # union divided by total edges
        inter = len([this for this in self.edges if 
            this in G.edges()]);
        return float(inter)/len(self.edges)

    def compute_distances(self, N=100, show=True, progress=True):
        self.init()
        from easydev import progress_bar
        distances = []
        pb = progress_bar(N)
        for i in range(0,N):
            self.swap(1, inplace=True)
            dist = self.get_distance(self.graph)
            distances.append(dist)
            if progress:pb.animate(i)
        if show is True:
            import pylab
            pylab.plot(distances)
            pylab.grid(True)
        return distances

    def plot_average_distance(self, repeat=10, N=100):
        import pandas as pd
        import pylab
        distances = []
        
        from easydev import progress_bar
        pb = progress_bar(repeat)
        for i in range(0, repeat):
            distance = self.compute_distances(N=N, show=False, progress=False)
            distances.append(distance)
            pb.animate(i+1)
        df = pd.DataFrame(distances)
        pylab.clf()
        pylab.fill_between(range(0,N), df.mean()+df.std(),
                y2=df.mean()-df.std()); 
        pylab.plot(df.mean(), 'r', lw=2)
        pylab.grid(True)
        pylab.ylim([0,1])
        pylab.ylabel('similarity')
        return distances
