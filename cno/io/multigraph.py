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

import pylab
import networkx as nx

# cellnopt modules
from cno.io.sif import SIF
from cno.io.midas import XMIDAS
from cno.io.cnograph import CNOGraph

from colormap import Colormap

__all__ = ["CNOGraphMultiEdges"]




class CNOGraphMultiEdges(CNOGraph, nx.MultiDiGraph):
    """A multiDigraph version of CNOGraph.

    .. warning:: experimental

    .. plot::
        :width: 50%
        :include-source:

         >>> from cno import cnodata
         >>> from cno.io.multigraph import CNOGraphMultiEdges
         >>> c = CNOGraphMultiEdges(cnodata("PKN-ToyPCB.sif"), cnodata("MD-ToyPCB.csv"))
         >>> c.add_edge("PAK", "Mek", link="-")
         >>> c.plot()

    .. plot::
        :width: 50%
        :include-source:

        >>> from cno.io.multigraph import CNOGraphMultiEdges
        >>> c2 = cnograph.CNOGraphMultiEdges()
        >>> c2.add_edge("A","B", link="+", edgecolor=.1)
        >>> c2.add_edge("A","B", link="+", edgecolor=.2)
        >>> c2.add_edge("A","B", link="-", edgecolor=.3)
        >>> c2.add_edge("A","B", link="-", edgecolor=.4)
        >>> c2.plot(edge_attribute="edgecolor", colorbar=True, cmap="spring")


    .. todo:: self.reacID attribute does not work
    """
    def __init__(self, model=None, data=None, verbose=False, **kargs):
        super(CNOGraphMultiEdges, self).__init__(**kargs)
        self._graph_type = 'multigraph'
        c = CNOGraph(model, data)
        self.midas = c.midas
        for e in c.edges(data=True):
            self.add_edge(e[0], e[1], **e[2])

    def reset_edge_attributes(self):
        """set all edge attributes to default attributes

        required to overwrite the cnograph method that do to handle the multigraph structure

        .. seealso:: :meth:`~cno.io.cnograph.CNOGraph.set_default_edge_attributes`

        """
        for edge in self.edges():
            for key in self.edge[edge[0]][edge[1]].keys():

                attrs = self.edge[edge[0]][edge[1]][key]
                attrs = self.set_default_edge_attributes(**attrs)
                self.edge[edge[0]][edge[1]][key] = attrs

    def _set_edge_attribute_color(self, edge_attribute, cmap):
        """Need to overwrite the cnograph method to handle the multigraph structure"""
        import matplotlib
        cmap = matplotlib.cm.get_cmap(cmap)
        sm = matplotlib.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(vmin=0,vmax=1), cmap=cmap)
        M = max([self.edge[edge[0]][edge[1]][key][edge_attribute] for edge in self.edges() for key in self.edge[edge[0]][edge[1]].keys()])
        for edge in self.edges():
            for key in self.edge[edge[0]][edge[1]].keys():
                value = self.edge[edge[0]][edge[1]][key][edge_attribute]/float(M)
                rgb = sm.to_rgba(value)
                colorHex = matplotlib.colors.rgb2hex(rgb)
                self.edge[edge[0]][edge[1]][key]['color'] = colorHex

    def _ambiguous_multiedge(self, node):
        """This is a multiedge so it is always ambiguous ?
        """
        return True

    def is_compressable(self, node):
        return False

    #def _get_compressable_nodes(self):
    #    return []
    #def _get_compressable(self):
    #    return []
    #compressable_nodes = property(_get_compressable)

    def remove_edge(self, u, v, key=None):
        """Remove the edge between u and v.

        :param str u: node u
        :param str u: node v

        Calls :meth:`clean_orphan_ands` afterwards
        """
        super(CNOGraph, self).remove_edge(u,v, key=key)
        #if "+" not in n:
        self.clean_orphan_ands()

    def _set_edge_attribute_label(self, this, edge_attribute):
        for e in this.edges():
            for key in self.edge[e[0]][e[1]].keys():
                this.edge[e[0]][e[1]][key]['label'] = this.edge[e[0]][e[1]][key][edge_attribute]







