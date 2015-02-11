# -*- python -*-
#
#  This file is part of the cinapps.tcell package
#
#  Copyright (c) 2012-2013 - EMBL-EBI
#
#  File author(s): Thomas Cokelaer (cokelaer@ebi.ac.uk)
#
#  Distributed under the GLPv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: www.cellnopt.org
#
##############################################################################
from __future__ import print_function

from matplotlib import colors
import pylab
import networkx as nx
import numpy as np
import pandas as pd

from cno.io.cnograph import CNOGraph

__all__ = ["XCNOGraph"]



class XCNOGraph(CNOGraph):
    """Extra plotting and statistical tools related to CNOGraph"""
    def __init__(self, model=None, midas=None, verbose=False):
        super(XCNOGraph, self).__init__(model, midas, verbose=verbose)

    def hcluster(self):
        """

        .. plot::
            :include-source:
            :width: 50%

            from cno import XCNOGraph, cnodata
            c = XCNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
            c.hcluster()

        .. warning:: experimental
        """
        from scipy.cluster import hierarchy
        from scipy.spatial import distance
        path_length=nx.all_pairs_shortest_path_length(self.to_undirected())
        n = len(self.nodes())
        distances=np.zeros((n,n))
        nodes = self.nodes()
        for u,p in path_length.items():
            for v,d in p.items():
                distances[nodes.index(u)-1][nodes.index(v)-1] = d
        sd = distance.squareform(distances)
        hier = hierarchy.average(sd)
        pylab.clf();
        hierarchy.dendrogram(hier)

        pylab.xticks(pylab.xticks()[0], nodes)

    def plot_degree_rank(self, loc='upper right', alpha=0.8, markersize=10,
            node_size=25, layout='spring', marker='o', color='b'):
        """Plot degree of all nodes

        .. plot::
            :include-source:
            :width: 50%

            from cno import XCNOGraph, cnodata
            c = XCNOGraph(cnodata("PKN-ToyPB.sif"))
            c.plot_degree_rank()

        """
        degree_sequence=sorted(nx.degree(self).values(),reverse=True) # degree sequence

        pylab.clf()
        pylab.loglog(degree_sequence, color+'-', marker=marker,
                markersize=markersize)
        pylab.title("Degree/rank and undirected graph layout")
        pylab.ylabel("Degree")
        pylab.xlabel("Rank")

        # draw graph in inset
        if loc == 'upper right':
            pylab.axes([0.45, 0.45, 0.45, 0.45])
        else:
            pylab.axes([0.1, 0.1, 0.45, 0.45])

        UG = self.to_undirected()
        Gcc = list(nx.connected_component_subgraphs(UG))
        Gcc = Gcc[0]
        if layout == 'spring':
            pos = nx.spring_layout(Gcc)
        else:
            pos = nx.circular_layout(Gcc)
        pylab.axis('off')
        nx.draw_networkx_nodes(Gcc, pos, node_size=node_size)
        nx.draw_networkx_edges(Gcc, pos, alpha=alpha)
        pylab.grid()
        #pylab.show()

    def plot_feedback_loops_histogram(self, **kargs):
        """Plots histogram of the cycle lengths found in the graph

        :return: list of lists containing all found cycles
        """
        data = list(nx.simple_cycles(self))
        if len(data):
            pylab.hist([len(x) for x in data], **kargs)
            pylab.title("Length of the feedback loops")
        else:
            print('No loop/cycle found')
        return data

    def plot_in_out_degrees(self, show=True,ax=None, kind='kde'):
        """
         .. plot::
            :include-source:
            :width: 50%

            from cno import XCNOGraph, cnodata
            c = XCNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
            c.plot_in_out_degrees()


        """
        ts1 = pd.TimeSeries(self.in_degree())
        ts2 = pd.TimeSeries(self.out_degree())
        df = pd.DataFrame([ts1, ts2]).transpose()
        df.columns = ["in","out"]
        if show:
            df.plot(kind=kind, ax=ax)  # kernerl density estimation (estimiation of histogram)
        #df = ...
        #df.transpose().hist()
        return df

    def plot_feedback_loops_species(self, cmap="heat", **kargs):
        """Returns and plots species part of feedback loops


        :param str cmap: a color map
        :return: dictionary with key (species) and values (number of feedback loop
            containing the species) pairs.

        .. plot::
            :include-source:
            :width: 50%

            from cno import XCNOGraph, cnodata
            c = XCNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
            c.plot_feedback_loops_species(cmap='heat', colorbar=True)

        """
        if len(self) == 0:
            self.logging.warning("Empty graph")
            return

        data = nx.simple_cycles(self)
        data = list(pylab.flatten(data))

        # FIXME: may not be robust to have "and": could be a valid name
        counting = [(x, data.count(x)) for x in self.nodes() 
                if data.count(x)!=0 and str(x).startswith('and') is False 
                and self.isand(x) is False]

        for node in self.nodes():
            self.node[node]['loops'] = 0

        if len(counting):
            M = float(max([count[1] for count in counting]))
            # set a default
            #for node in self.nodes():
            #    self.node[node]['loops'] = "#FFFFFF"

            for count in counting:
                #ratio_count = sm.to_rgba(count[1]/M)
                ratio_count = count[1]/M
                colorHex = ratio_count
                #self.node[count[0]]['loops'] = colorHex
                self.node[count[0]]['loops'] = ratio_count
                self.node[count[0]]['style'] =  'filled,bold'

        self.plot(node_attribute="loops", cmap=cmap, **kargs)
        return counting

    def degree_histogram(self, show=True, normed=False):
        """Compute histogram of the node degree (and plots the histogram)

        .. plot::
            :include-source:
            :width: 50%

            from cno import XCNOGraph, cnodata
            c = XCNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
            c.degree_histogram()


        """
        degree = self.degree().values()
        Mdegree = max(degree)

        if show == True:
            pylab.clf()
            res = pylab.hist(degree, bins=range(0,Mdegree+1), align='left',
                             rwidth=0.8, normed=normed)
            xlims = pylab.xlim()
            ylims = pylab.ylim()
            pylab.axis([0, xlims[1], ylims[0], ylims[1]*1.1])
            pylab.grid()
            pylab.title("Degree distribution")
        return res

    def plot_adjacency_matrix(self, fontsize=12, **kargs):
        """Plots adjacency matrix

        :param kargs: optional arguments accepted by pylab.pcolor

        From the following graph, 

        .. plot::
            :width: 70%

            from cno import XCNOGraph, cnodata
            c = XCNOGraph(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
            c.plot()

        The adjacency matrix can be created as follows:

        .. plot::
            :width: 70%
            :include-source:

            from cno import XCNOGraph, cnodata
            c = XCNOGraph(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
            c.plot_adjacency_matrix()

        """
        nodeNames = sorted(self.nodes())
        nodeNamesY = sorted(self.nodes())

        nodeNamesY.reverse()
        N = len(nodeNames)

        data = self.adjacency_matrix(nodelist=nodeNames)
        # This is now a sparse matrix (networkx 1.9).
        try:
            data = data.todense()
        except:
            pass

        pylab.pcolor(pylab.flipud(pylab.array(data)), edgecolors="k", **kargs)
        pylab.axis([0, N, 0, N])
        pylab.xticks([0.5+x for x in pylab.arange(N)], nodeNames, rotation=90,
                      fontsize=fontsize)
        pylab.yticks([0.5+x for x in pylab.arange(N)], nodeNamesY, rotation=0,
                      fontsize=fontsize)
        try:pylab.tight_layout()
        except:pass

    def dependency_matrix(self, fontsize=12):
        r"""Return dependency matrix

        * :math:`D_{i,j}` = green ; species i is an activator of species j (only positive path)
        * :math:`D_{i,j}` = red   ; species i is an inhibitor of species j (only negative path)
        * :math:`D_{i,j}` = yellow; ambivalent (positive and negative paths connecting i and j)
        * :math:`D_{i,j}` = red   ; species i has no influence on j

        .. plot::
            :include-source:
            :width: 80%

            from cno import XCNOGraph, cnodata
            c = XCNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
            c.dependency_matrix()

        """
        nodes = sorted(self.nodes())
        N = len(nodes)
        data = np.zeros((len(nodes), len(nodes)))
        for i,node1 in enumerate(nodes):
            paths = nx.shortest_path(self, node1)
            for j,node2 in enumerate(nodes):
                if node1 == node2:
                    data[i][j] = 0
                elif node2 not in paths.keys():
                    data[i][j] = 0
                else:
                    path = paths[node2]
                    links = [self.edge[path[ii]][path[ii+1]]["link"] for ii in range(0,len(path)-1)]
                    if len(np.unique(links)) == 2:
                        data[i][j] = 1  # yellow
                    elif "+" in links:
                        data[i][j] = 2  #green
                    elif "-" in links:
                        if links.count("-") % 2 ==0:
                            data[i][j] = 2
                        else:
                            data[i][j] = 3   #red

        nodeNames = [node.replace("_", "\_") for node in nodes]
        nodeNamesY = [node.replace("_", "\_") for node in nodes]

        norm = colors.Normalize(vmin=0, vmax=3)

        cmap = colors.ListedColormap([[0., 0, 0], [1,1,0],[.5, 1, 0.], [1., 0, 0.]])
        indices = [i for i, node in enumerate(nodes) 
                if "and" not in node or "+" in nodes]

        pylab.clf()
        pylab.pcolor(pylab.flipud(data[indices][:,indices]), edgecolors="w", 
                cmap=cmap, norm=norm);

        N = len(indices)

        nodeNames = np.array(nodeNames)[indices]
        nodeNamesY = np.array(nodeNamesY)[indices[::-1]]
        X = [0.5+x for x in range(0, len(nodeNames))]
        pylab.xticks(X,  nodeNames, rotation=90, fontsize=fontsize)
        pylab.yticks(X,  nodeNamesY, fontsize=fontsize)
        pylab.xlim([0, len(X)])
        pylab.ylim([0, len(X)])


    def random_cnograph(self, nStimuli=5, nSignals=14, fraction_activation=0.9, nTranscript=5, 
            nExtraNode=10):
        """

        Create the nodes first (stimuli, signals, transcripts, extra nodes). Them
        add edges such that the ratio of activation/inhibition is fixed.

        no self loop
        """
        assert fraction_activation >=0 and fraction_activation<=1
        assert nStimuli>=1
        assert nExtraNode >= 1
        assert nSignals >= 1
        assert nTranscript >=1 and nTranscript <= nSignals

        self.clear()

        # add stimuli
        stimuli = ['L' + str(i) for i in range(1, nStimuli+1)]
        self.add_nodes_from(stimuli)
        self._stimuli = stimuli[:]

        signals = ['S' + str(i) for i in range(1, nSignals+1)]
        self.add_nodes_from(signals)
        self._signals = signals[:]
        
        self.add_nodes_from(['N'+str(i) for i in range(1,nExtraNode+1)])
        nodes = self.nodes()

        # select the transcript:
        transcripts = [x for x in self.nodes() if x not in self.stimuli]
        transcripts = transcripts[0:nTranscript]

        def get_link():
            link = np.random.uniform()
            if link < fraction_activation:
                return "+"
            else:
                return "-"

        count = 0
        N = len(self.nodes())
        while nx.is_connected(self.to_undirected()) is False and count < N * 3:
            np.random.shuffle(nodes)
            n1 = nodes[0]
            n2 = nodes[1]
            if n2 in self.stimuli or n1 in transcripts:
                continue
            # ignore self loop
            if n1 == n2:
                continue
            self.add_edge(n1, n2, link=get_link())
            count += 1


        # some nodes (non-stimuli/non-signals) may have no input connections, which is
        # not wanted.
        tofix = [x for x in self.nodes() if self.in_degree()[x] == 0 and x.startswith('N')]
        for nodetofix in tofix:
            nodes = [node for node in self.nodes() if node !=tofix]
            np.random.shuffle(nodes)
            self.add_edge(nodes[0], nodetofix, link=get_link())

        # make sure the ligands are connected:
        for stimulus in self._stimuli:
            if len(self.successors(stimulus)) == 0:
                print("fixing stimulus %s" % stimulus)
                nodes = [node for node in self.nodes() if node not in self._stimuli]
                np.random.shuffle(nodes)
                self.add_edge(stimulus, nodes[0])


    def random_poisson_graph(self, n=10, mu=2.5, ratio=0.9, 
            remove_unconnected=True, Nsignals=5, Nstimuli=5, 
            remove_self_loops=True, maxtrials=50):
        """Experimental random graph creation"""
        count = 0
        while count < maxtrials:
            self._random_poisson_graph(n, mu, ratio=ratio,
                remove_unconnected=remove_unconnected, 
                remove_self_loops=remove_self_loops)
            if nx.is_connected(self.to_undirected()):
                count = maxtrials + 1
            else:
                count += 1

    def _random_poisson_graph(self, n=10, mu=2.5, ratio=0.9, 
            remove_unconnected=True, 
            remove_self_loops=True,  Nsignals=5, Nstimuli=5):
        from scipy.stats import poisson
        z = [poisson.rvs(mu) for i in range(0,n)]
        G = nx.expected_degree_graph(z)
        self.clear()

        # converts to strings
        edges = [(str(e[0]), str(e[1])) for e in G.edges()]
        assert ratio >= 0
        assert ratio <= 1

        N = int(len(edges)* ratio)
        edges_pos = edges[0:N]
        edges_neg = edges[N:]
        self.add_edges_from(edges_pos, link="+")
        self.add_edges_from(edges_neg, link="-")

        # remove self loop first
        if remove_self_loops:
            self.remove_self_loops()

        if remove_unconnected == False:
            # add all nodes (even though they me be unconnected
            self.add_nodes_from(G.nodes())

        ranks = self.get_same_rank()
        sources = ranks[0]
        sinks = ranks[max(ranks.keys())]
        Nstim = min(len(sources), Nstimuli)
        Nsignals = min(len(sinks), Nsignals)
        self._stimuli = sources[0:Nstim]
        self._signals = sinks[0:Nsignals]
        self.set_default_node_attributes()
