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
""".. topic:: **One of the main data structures of cellnopt to manipulate networks**"""
from __future__ import print_function
import os
import copy
import tempfile
import itertools
import subprocess
import shutil
import json
import collections
from functools import wraps

import matplotlib
import pylab
import networkx as nx

from cno.misc.profiler import do_profile

try:
    import pygraphviz as gv
except ImportError:
    print("Warning:: Pygraphviz not found")
    pass

import numpy as np
from easydev import Logging, TempFile, DevTools

# cellnopt modules
from cno.io.sif import SIF
from cno.io.midas import XMIDAS
from cno.io.reactions import Reaction
from cno.misc import CNOError
import colormap

from cno.misc.profiler import do_profile
__all__ = ["CNOGraph", "CNOGraphAttributes"]



def modifier(func):
    @wraps(func)
    def wrapper(self, *args, **kargs):
        return func(self, *args, **kargs)
    return wrapper


class Attributes(dict):
    """Simple dictionary to handle attributes (nodes or eges)"""
    def __init__(self, color, **kargs):
        self['color'] = color
        for k,v in kargs.items():
            self[k] = v


class EdgeAttributes(Attributes):
    """Simple dictionary to handle edge attributes"""
    def __init__(self, penwidth=1, color="black", arrowhead="normal", **kargs):
        super(EdgeAttributes, self).__init__(color=color, **kargs)
        self['penwidth'] = penwidth
        self['arrowhead'] = arrowhead


class NodeAttributes(dict):
    """Simple dictionary to handle node attributes

    Used by the :class:`CNOGraph`.

    ::

        >>> a = Attributes(color="red", fillcolor="white")
        >>> assert a['color'] == "red"
        True

    """
    def __init__(self, color="black", fillcolor="white", shape="rectangle",
                 style="filled,bold", **kargs):
        """

        :param color: color of the border
        :param fillcolor: color inside the shape
        :param shape color inside the shape
        :param style: color inside the shape
        """
        super(NodeAttributes, self).__init__(color=color, **kargs)

        #self['color'] = color
        self['fillcolor'] = fillcolor
        self['shape'] = shape
        self['style'] = style


class CNOGraphAttributes(object):
    """

    define attributes of the cnograp when calling the plotting function.
    """
    def __init__(self):

        self.attributes = {
               'activation': EdgeAttributes(color='black', arrowhead="normal"),
               'inhibition': EdgeAttributes(color='red', arrowhead="tee"),
               'stimuli': NodeAttributes(fillcolor='#9ACD32'),
               'inhibitors': NodeAttributes(fillcolor='orangered'),
               'signals': NodeAttributes(fillcolor='lightblue'),
               'nonc': NodeAttributes( fillcolor='gray', style='diagonals, filled'),
               'compressed': NodeAttributes(fillcolor="white",style='dashed', penwidth=2),
               'others': NodeAttributes(fillcolor="white",style='filled,bold',penwidth=2),
               'and': NodeAttributes(shape='circle',style='filled',width=.1,
                            height=.1, fixedsize=True, label='')
               }
    def keys(self):
        return self.attributes.keys()
    def __getitem__(self, key):
        return self.attributes[key]


class CNOGraph(nx.DiGraph):
    """Data structure (Digraph) used to manipulate networks

    The networks can represent for instance a protein interaction network (PIN).

    CNOGraph is a data structure dedicated to the analysis of
    phosphorylation data within protein-protein interaction networks
    but can be used in a more general context. Note that CNOGraph inherits
    from the **directed graph** data structure of networkx.

    However, we impose links between nodes to be restricted to two types:
        * "+" for activation
        * "-" for inhibition.

    An empty instance can be created as follows::

        c = CNOGraph()

    and edge can be added as follows::

        c.add_edge("A", "B", link="+")
        c.add_edge("A", "C", link="-")

    An even simpler way is to add :class:`cno.io.reactions.Reaction`, which can be strings
    or instance of the Reaction class.

    The methods :meth:`add_node`  and :meth:`add_edge` methods can be used to
    populate the graph. However, it is also possible to read a network
    stored in a file in :class:`cno.io.sif.SIF` format::

        >>> from cno import CNOGraph, cnodata
        >>> pknmodel = cnodata("PKN-ToyPB.sif")
        >>> c = CNOGraph(pknmodel)

    The SIF model can be a filename, or an instance of
    :class:`~cno.io.sif.SIF`. Note for CellNOpt users
    that if **and** nodes are contained in the original SIF files, they are
    transformed int AND gates using "^" as the logical AND.

    Other imports are available, in particular :meth:`read_sbmlqual`.

    You can add or remove nodes/edges in the CNOGraph afterwards using NetworkX methods.

    When instanciating a CNOGraph instance, you can also populate data from
    a :class:`~cno.io.midas.XMIDAS` data instance or a MIDAS filename.
    MIDAS file contains measurements made on proteins
    in various experimental conditions (stimuli and inhibitors). The names of
    the simuli, inhibitors and signals are used to color the nodes in the
    plotting function. However, the data itself is not used.

    If you don't use any MIDAS file as input, you can set the
    stimuli/inhibitors/signals manually by filling the hidden attributes
    _stimuli, _signals and _inhibitors with list of nodes contained in the graph.

    .. rubric:: Node and Edge attributes

    The node and edge attributes can be accessed as follows (and changed)::

        >>> c.node['egf']
        {'color': u'black',
         u'fillcolor': u'white',
         'penwidth': 2,
         u'shape': u'rectangle',
         u'style': u'filled,bold'}

        >>> c.edge['egf']['egfr']
        {u'arrowhead': u'normal',
         u'color': u'black',
         u'compressed': [],
         'link': u'+',
         u'penwidth': 1}

    .. rubric:: OPERATORS

    CNOGraph is a data structure with useful operators (e.g. union). Note,
    however, that these operators are applied on the topology only (MIDAS
    information is ignored). For instance, you can add graphs with the **+** operator
    or check that there are identical ::

        c1 = CNOGraph()
        c1.add_reaction("A=B")
        c2 = CNOGraph()
        c2.add_reaction("A=C")
        c3 = c1 +c2

    Let us illustrate the + operation with another example. Let us consider the following graphs:

    .. plot::
        :include-source:
        :width: 30%

        from cno import CNOGraph
        c1 = CNOGraph()
        c1.add_edge("A","B", link="+")
        c1.add_edge("A","C", link="-")
        c1.plot()


    .. plot::
        :include-source:
        :width: 30%

        from cno import CNOGraph
        c2 = CNOGraph()
        c2.add_edge("A","E", link="+")
        c2.add_edge("C","E", link="+")
        c2.plot()

    ::

        (c1+c2).plot()


    .. plot::
        :width: 50%

        from cno import CNOGraph
        c1 = CNOGraph()
        c1.add_edge("A","B", link="+")
        c1.add_edge("A","C", link="-")
        c1.plot()
        c2 = CNOGraph()
        c2.add_edge("A","E", link="+")
        c2.add_edge("C","E", link="+")
        c2.plot()
        (c1+c2).plot()

    You can also substract a graph from another one::

        c3 = c1 - c2
        c3.nodes()

    The new graph should contains only one node (B). Additional functionalities
    such as :meth:`intersect`, :meth:`union` and :meth:`difference` can be used to see the difference
    between two graphs.

    .. rubric:: PLOTTING

    There are plotting functionalities to look at the graph, which are based on graphviz
    library. For instance, the :meth:`plot` function is quite flexible. If a MIDAS file
    is provided, the default behaviour follow CellNOptR convention,  where stimuli are
    colored in green, inhibitors in red and measurements in blue:

    .. plot::
        :include-source:
        :width: 50%

        from cno import CNOGraph, cnodata
        pknmodel = cnodata("PKN-ToyPB.sif")
        data = cnodata("MD-ToyPB.csv")
        c = CNOGraph(pknmodel, data)
        c.plot()

    If you did not use any MIDAS file as input parameter, you can still populate the hidden fields
    :attr:`_stimuli`, :attr:`_inhibitors`, :attr:`_signals`.

    You can also overwrite this behaviour by using the node_attribute parameter when
    calling :meth:`plot`. For instance, if you call :meth:`centrality_degree`, which
    computes and populate the node attribute
    **degree**. You can then call plot as follows to replace the default
    color:

    .. plot::
        :include-source:
        :width: 50%

        from cno import CNOGraph, cnodata
        pknmodel = cnodata("PKN-ToyPB.sif")
        data = cnodata("MD-ToyPB.csv")
        c = CNOGraph(pknmodel, data)
        c.centrality_degree()
        c.plot(node_attribute="centrality_degree", colorbar)

    Similarly, you can tune the color of the edge attribute. See the :meth:`plot` for more details.

    .. seealso::  tutorial, user guide

    .. seealso:: The :class:`cno.io.xcnograph.XCNOGraph` provides many more tools for plotting
        various information on the graph structure.


    """
    def __init__(self, model=None, data=None, verbose=False, **kargs):
        """.. rubric:: Constructor

        :param str model: optional network in SIF format. Can be the filename
            or instance of :class:`~cno.io.sif.SIF`
        :param data: optional data file in MIDAS format. Can be a filename or
            instance of :class:`~cno.io.midas.XMIDAS`
        :param bool verbose:
        :param str celltype: if a MIDAS file contains more that 1 celltype, you
            must provide a celltype name


        """
        super(CNOGraph, self).__init__(**kargs)
        self.kargs = kargs.copy()

        self._graph_type = 'digraph'
        # This is a DIgraph attribute
        # self.graph is a DiGraph attribute that is overwritten sometinmes

        self.graph_options = {
           'graph': {
                "title": "CNOGraph output generated by cno",
                "dpi":200,
                'rankdir':'TB', # TB, LR, RL
                #'ordering': "out",
                'splines':True,
                'fontsize': 22,
                 #'nodesep': .5,
                 #'ranksep':.6,
                'ratio':'auto', # numeric,  'fill','compress','expand','auto'
                # 0.8 is good for laptop screens. 2 is good for
                'size': "15,15",
                # 'fontname': 'helvetica',
                },
            'node':{
                #'width':1,
                #'fontsize':40,
                #'height':1,
                #'width':2,
                'fontname':'bold'
                 },
            'edge': {
                'minlen':1,
                'color':'black'
                },
            'ipython': {
                'width': 500,
                'rank_method': 'cno'
                }
            }

        # cellnoptR has always the same layout:
        #s.model.graph_options['graph']['nodesep'] = 0.5
        #s.model.plot(rank_method='same')

        self.plot_options = {
                'colorbar.orientation': 'horizontal',
                'colorbar.shrink': 0.5,
                'colorbar.fraction': 0.15,
                'colorbar.pad': 0.1,
                }

        #: nodes and edges attributes. See :class:`CNOGraphAttributes`
        self.attributes = CNOGraphAttributes()
        self._set_default_node_attribute = 'cno'


        self.and_symbol = "^"
        self.or_symbol = "+"

        self._midas = None
        self._verbose = verbose
        self.logging = Logging(self.verbose)

        #: stimuli
        self._stimuli = []
        #: inhibitors
        self._inhibitors = []

        self._compressed = []
        self._signals =[]
        self._nonc = None
        self._ranks = None

        # the model
        if hasattr(model, '__class__') and \
            model.__class__.__name__ in ['CNOGraph', 'XCNOGraph']:
            for node in model.nodes():
                self.add_node(str(node))
            for edge in model.edges(data=True):
                if "link" in edge[2]:
                    self.add_edge(str(edge[0]), str(edge[1]), link=edge[2]['link'])
                else:
                    self.add_edge(str(edge[0]), str(edge[1]), link="+")
            self.set_default_node_attributes() # must be call if sif or midas modified.
            self.filename = None
            if model.midas is not None:
                self.midas = model.midas.copy()
        elif model is None:
            self.filename = None
        elif isinstance(model, str):
            if model.endswith('.sif'):
                self.read_sif(model)
                self.filename = model[:]
            elif model.endswith(".xml"):
                self.read_sbmlqual(model)
                self.filename = model[:]
            else:
                raise CNOError("Only filenames with .sif and .xml (SBML-qual) extension are recognised.")
        elif isinstance(model, SIF):
            self.read_sif(model)
            self.filename = 'undefined'

        # the data
        if self.midas is None:
            self.midas = data

        self._colormap = colormap.Colormap()
        self._changed = True

    def _set_verbose(self, verbose):
        self.logging.debugLevel = verbose
        self.midas.logging.debugLevel = verbose
    def _get_verbose(self):
        return self._verbose
    verbose = property(_get_verbose, _set_verbose)

    # SOME PROPERTIES
    def _get_midas(self):
        return self._midas
    def _set_midas(self, data):
        self._changed = True
        if isinstance(data, str):
            self._midas = XMIDAS(data, cellLine=self.kargs.get("cellLine", None),
                                 verbose=self.verbose)
        elif isinstance(data, XMIDAS):
            self._midas = copy.deepcopy(data)
        elif data is None:
            self._midas = data
        else:
            msg = "Incorrect data, Must a valid MIDAS file or instance of XMIDAS class {}"
            raise ValueError(msg.format(data))
        self.check_data_compatibility()
        self.set_default_node_attributes()
    midas = property(fget=_get_midas, fset=_set_midas,
                     doc="MIDAS Read/Write attribute.")

    def _get_stimuli(self):
        stimuli = list(self._stimuli[:])
        if self.midas:
            stimuli += self.midas.names_stimuli[:]
        return stimuli
    stimuli = property(_get_stimuli,
            doc="list of stimuli found in the :attr:`midas` and hidden attribute :meth:`_stimuli`")

    def _get_inhibitors(self):
        inhibitors = list(self._inhibitors)
        if self.midas:
            inhibitors += self.midas.names_inhibitors[:]
        return inhibitors
    inhibitors = property(_get_inhibitors,
            doc="list of inhibitors found in the :attr:`midas` and hidden attribute :attr:`_inhibitors`")

    def _get_signals(self):
        signals = list(self._signals) # do not reference
        if self.midas:
            signals += list(self.midas.names_signals)
        return signals
    signals = property(_get_signals,
            doc="list of signals found in the :attr:`midas` and hidden attribute :meth:`_signals`")

    # METHODS
    def read_sif(self, model):
        """If the SIF file changes, we need to rebuild the graph."""
        # takes the SIF input file and build up the CNOGraph. remove all nodes
        # before
        self.clear()
        self.logging.debug("Reading the model")
        #if isinstance(model, (str, unicode)):
        if isinstance(model, (str)):
            sif = SIF(model)
        elif isinstance(model, SIF):
            sif = model
        else:
            raise ValueError("The sif input must be a filename to a SIF file or an instance of the SIF class")

        # add all reactions
        self.add_reactions(sif.reactions)

        # now, we need to set the attributes, only if we have a cnolist,
        # otherwise color is the default (white)
        self.set_default_node_attributes() # must be call if sif or midas modified.
        self.logging.debug("Model loaded")

    def _add_simple_reaction(self, reac, node_dict=None, edge_dict=None):
        """A=B or !A=B"""

        #reac = Reaction(reac) # validate the reaction
        #reac = reac.name
        self._changed = True
        this = self.attributes.attributes['others']
        if node_dict:
            this.update(node_dict)
        node_dict = this

        lhs, rhs = reac.split("=", 1)
        if rhs == "":
            self.add_node(lhs, attr_dict=node_dict)
        else:
            if lhs.startswith("!"):
                link = "-"
                lhs = lhs[1:]
            else:
                link = "+"
            if self.has_edge(lhs, rhs):
                if self[lhs][rhs]['link'] == link:
                    self.logging.info("skip existing reactions %s %s %s" % (lhs, link, rhs))
                else:
                    self.add_edge(lhs, rhs, attr_dict=edge_dict, link=link)
            else:
                self.add_edge(lhs,rhs, attr_dict=edge_dict, link=link)
        self._changed = True

    def add_reactions(self, reactions, node_dict=None, edge_dict=None):
        self._changed = True
        for reac in reactions:
            self.add_reaction(reac, node_dict=node_dict, edge_dict=edge_dict)

    def is_reaction_in(self, reaction):
        """Return boolean with present of a reaction in the network

        :param reaction: a valid reaction

        .. todo:: Could be a string of reaction

        """
        reaction = Reaction(reaction)
        if reaction in self.reactions:
            return True
        else:
            return False

    @modifier
    def add_reaction(self, reac, node_dict=None, edge_dict=None):
        """Add nodes and edges given a reaction

        :param str reac: a valid reaction. See below for examples

        Here are some valid reactions that includes NOT, AND and OR gates. + is an OR
        and ^ character is an AND gate::

            >>> s.add_reaction("A=B")
            >>> s.add_reaction("A+B=C")
            >>> s.add_reaction("A^C=E")
            >>> s.add_reaction("!F+G=H")

        .. plot::
            :width: 50%
            :include-source:

            from cno import CNOGraph
            c = CNOGraph()
            c.add_reaction("a+b^c+e+d^h=Z")
            c.plot()

        .. warning:: component of AND gates are ordered alphabetically.
        """
        # add the nodes first so that the attributes are ste properly
        # make sure that (1) there is no extra spaces and (2) this is a string, not
        # a unicode
        reac = Reaction(str(reac.strip()))
        # no need to sort the reaction here since we will split the ORs and sort the 
        # AND gates.

        # make sure that nodes are created with the default (or use attributes)
        # before creating the edges.
        for node in reac.species:
            self.add_node(node, attr_dict=node_dict)

        # if there is an OR gate, easy, just need to add simple reactions
        # A+!B=C is split into A=C and !B=C

        for this_lhs in reac.lhs.split("+"):
            # + has priority upon ^ unlike in maths so we can split with +
            # A+B^C^D+E=C means 3 reactions: A=C, E=C and B^C^D=C
            if self.isand(this_lhs) is False:
                name = this_lhs + "=" + reac.rhs
                self._add_simple_reaction(name, node_dict=node_dict,
                        edge_dict=edge_dict)
            else:
                and_gate_name = this_lhs + "=" + reac.rhs
                reac = Reaction(and_gate_name)
                reac.sort()
                # and gates need a little bit more work
                # the and gate to the output
                self.add_edge(reac.name, reac.rhs, attr_dict=edge_dict,
                        link="+")
                # now the inputs to the and gate
                for this in this_lhs.split(self.and_symbol):
                    name = this + "=" + reac.name
                    self._add_simple_reaction(name,
                            node_dict=node_dict, edge_dict=edge_dict)
        self._changed = True

    #@do_profile()
    def get_default_edge_attributes(self,  **attr):
        #if "compressed" not in attr.keys():
        #    attr["compressed"] = []

        link = attr.get("link", "+")
        if link == '+':
            attr['color'] = 'black'
            attr['arrowhead'] = 'normal'
            attr['penwidth'] = attr.get('penwidth', 1)
        else:
            attr['color'] = 'red'
            attr['arrowhead'] = 'tee'
            attr['penwidth'] = attr.get('penwidth', 1)

        return attr

    def reset_edge_attributes(self):
        """set all edge attributes to default attributes

        .. seealso:: :meth:`get_default_edge_attribute`

        if we set an edge label, which is an AND ^, then plot fails in this function
        c.edge["alpha^NaCl=HOG1"]['label'] = "?"
        """
        for edge in self.edges():
            attrs = self.edge[edge[0]][edge[1]]
            attrs = self.get_default_edge_attributes(**attrs)
            self.edge[edge[0]][edge[1]] = attrs

    def add_edges_from(self, ebunch, attr_dict=None, **attr):
        self._changed = True
        for edge in ebunch:
            self.add_edge(edge[0], edge[1], attr_dict=attr_dict, **attr)

    #@do_profile()
    def add_edge(self, u, v, attr_dict=None, **attr):
        """adds an edge between node u and v.

        :param str u: source node
        :param str v: target node
        :param str link: compulsary keyword. must be "+" or "-"
        :param dict attr_dict: dictionary, optional (default= no attributes)
             Dictionary of edge attributes.  Key/value pairs will update existing
             data associated with the edge.
        :param attr: keyword arguments, optional
            edge data (or labels or objects) can be assigned using keyword arguments.
            keywords provided will overwrite keys provided in the **attr_dict** parameter

        .. warning:: color, penwidth, arrowhead keywords are populated according to the
            value of the link.

        * If link="+", then edge is black and arrowhead is normal.
        * If link="-", then edge is red and arrowhead is a tee

        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph
            c = CNOGraph()
            c.add_edge("A","B",link="+")
            c.add_edge("A","C",link="-")
            c.add_edge("C","D",link="+", mycolor="blue")
            c.add_edge("C","E",link="+", data=[1,2,3])

        If you want multiple edges, use add_reaction() method.

            c.add_reaction("A+B+C=D")

        equivalent to ::

            c.add_edge("A", "D", link="+")
            c.add_edge("B", "D", link="+")
            c.add_edge("C", "D", link="+")


        Attributes on the edges can be provided using the parameters **attr_dict** (a dictionary)
        and/or ****attr**, which is a list of key/value pairs. The latter will overwrite the
        key/value pairs contained in the dictionary. Consider this example::

            c = CNOGraph()
            c.add_edge("a", "c", attr_dict={"k":1, "data":[0,1,2]}, link="+", k=3)
            c.edges(data=True)
            [('a',
            'c',
                {'arrowhead': 'normal',
                'color': 'black',
                'compressed': [],
                'k':3
                'link': '+',
                'penwidth': 1})]

        The field "k" in the dictionary (attr_dict) is set to 1. However, it is also
        provided as an argument but with the value 3. The latter is the one used to populate
        the edge attributes, which can be checked by printing the data of the edge (c.edges(data=True())


        .. seealso:: special attributes are automatically set by :meth:`get_default_edge_attributes`.
            the color of the edge is black if link is set to "+" and red otherwie.

        """

        # cast u to str to search for + sign
        if "+" in str(u):
            lhs = u.split("+")
            for x in lhs:
                if x.startswith("!"):
                    attr["link"] = "-"
                    attr['color'] = 'red'
                    attr['arrowhead'] = 'tee'
                    super(CNOGraph, self).add_edge(x[1:], v, attr_dict, **attr)
                else:
                    attr["link"] = "+"
                    super(CNOGraph, self).add_edge(x, v, attr_dict, **attr)
        else:
            attr = self.get_default_edge_attributes(**attr)
            super(CNOGraph, self).add_edge(u, v, attr_dict, **attr)

    def clear(self):
        """Remove nodes and edges and MIDAS instance"""
        super(CNOGraph, self).clear()
        self._changed = True
        self.midas = None
        self._stimuli = []
        self._signals = []
        self._inhibitors = []
        self._ranks = None

    @modifier
    def clean_orphan_ands(self):
        """Remove AND gates that are not AND gates anymore

        When removing an edge or a node, AND gates may not be valid anymore
        either because the output does not exists or there is a single input.

        This function is called when :meth:`remove_node` or :meth:`remove_edge` are called.
        However, if you manipulate the nodes/edges manually you may need to call
        this function afterwards.
        """
        self._changed = True
        for node in self._find_and_nodes():
            if len(self.successors(node))==0 or len(self.predecessors(node))<=1:
                self.remove_node(node)
                continue

    def check_data_compatibility(self):
        """When setting a MIDAS file, need to check that it is compatible with
        the graph, i.e. species are found in the model."""
        if self.midas:
            msg = "The %s %s was found in the MIDAS file but is "
            msg += "not present in the model. Adding that node to your model"
            for x in self.midas.names_cues:
                if x not in self.nodes():
                    self.add_node(x)
                    self.logging.warning(msg % ('cues', x))
            for x in self.midas.names_signals:
                if x not in self.nodes():
                    self.add_node(x)
                self.logging.warning(msg % ('signals', x))

    @modifier
    def remove_and_gates(self):
        """Remove the AND nodes added by :meth:`expand_and_gates`"""
        self._changed = True
        for n in self._find_and_nodes():
            self.remove_node(n)

    def __eq__(self, other):
        # we must look at the data to figure out the link + or - but should ignore
        # all other keys
        edges1 = sorted(self.edges(data=True))
        edges2 = sorted(other.edges(data=True))
        edges1  = [(e[0], e[1], {'link':e[2]['link']}) for e in edges1]
        edges2  = [(e[0], e[1], {'link':e[2]['link']}) for e in edges2]
        res = edges1 == edges2
        return res

    def __add__(self, other):
        """allows a+b operation

        combines the _inhibitors, _signals, _stimuli but keep only the first
        midas file !

        """
        G = self.copy()
        G.add_nodes_from(other.nodes(data=True))
        edges = other.edges(data=True)
        for e1,e2,d in edges:
            G.add_edge(e1,e2,None, **d)
        G._inhibitors += other._inhibitors
        G._signals += other._signals
        G._stimuli += other._stimuli

        # TODO: merge the MIDAS files. ?
        return G

    def __sub__(self, other):
        G = self.copy()
        G.remove_nodes_from([n for n in G if n in other.nodes()])
        return G

    def __str__(self):
        nodes = len([x for x in self.nodes() if self.and_symbol not in str(x)])
        andnodes = len([x for x in self.nodes() if self.and_symbol in str(x)])

        msg = "The model contains %s nodes (and %s AND node)\n" % (nodes, andnodes)

        self.logging.warning("Edge counting valid only if and node have only 2 inputs")
        edges = len([e for e in self.edges() if self.and_symbol not in
            str(e[0]) and self.and_symbol not in str(e[1])])
        andedges = len([e for e in self.edges() if self.and_symbol
            in str(e[0]) or self.and_symbol in str(e[1])])/3
        msg += "%s Hyperedges found (%s+%s) \n" % (edges+andedges, edges, andedges)

        return msg

    #def __rsub__(self, other):
    #    self.remove_nodes_from([n for n in self if n in other.nodes()])

    def union(self, other):
        """Return graph with elements from this instance and the input graph.

        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph
            from pylab import subplot, title

            c1 = CNOGraph()
            c1.add_edge("A", "B", link="+")
            c1.add_edge("A", "C", link="-")
            c1.add_edge("C", "E", link="+")
            subplot(1,3,1)
            title(r"graph $C_1$")
            c1.plot(hold=True)

            c2 = CNOGraph()
            c2.add_edge("A", "B", link="+")
            c2.add_edge("B", "D", link="+")
            c2.add_edge("B", "F", link="+")
            subplot(1,3,2)
            c2.plot(hold=True)
            title(r"graph $C_2$")

            c3 = c1.union(c2)
            subplot(1,3,3)
            c3.plot(hold=True)
            title(r"graph $C_3 = C_1 \cup C_2$")

        """
        c = self + other
        return c

    def difference(self, other):
        """Return a CNOGraph instance that is the difference with the input graph

        (i.e. all elements that are in this set but not the others.)

        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph
            from pylab import subplot, title

            c1 = CNOGraph()
            c1.add_edge("A", "B", link="+")
            c1.add_edge("A", "C", link="-")
            c1.add_edge("C", "E", link="+")
            subplot(1,3,1)
            title("graph C1")
            c1.plot(hold=True)

            c2 = CNOGraph()
            c2.add_edge("A", "B", link="+")
            c2.add_edge("B", "D", link="+")
            c2.add_edge("B", "F", link="+")
            subplot(1,3,2)
            c2.plot(hold=True)
            title("graph C2")

            c3 = c1.difference(c2)
            subplot(1,3,3)
            c3.plot(hold=True)
            title("graph C3=C1-C2")


        .. note:: this method should be equivalent to the - operator. So c1-c2 == c1.difference(c2)
        """
        G = self.copy()
        G.remove_nodes_from([n for n in G if n in other.nodes()])
        return G

    def intersect(self, other):
        """Return a graph with only nodes found in "other" graph.

        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph
            from pylab import subplot, title

            c1 = CNOGraph()
            c1.add_edge("A", "B", link="+")
            c1.add_edge("A", "C", link="-")
            c1.add_edge("C", "E", link="+")
            subplot(1,3,1)
            title(r"graph $C_1$")
            c1.plot(hold=True)

            c2 = CNOGraph()
            c2.add_edge("A", "B", link="+")
            c2.add_edge("B", "D", link="+")
            c2.add_edge("B", "F", link="+")
            subplot(1,3,2)
            c2.plot(hold=True)
            title(r"graph $C_2$")

            c3 = c1.intersect(c2)
            subplot(1,3,3)
            c3.plot(hold=True)
            title(r"graph $C_3 = C_1 \cap C_2$")

        """
        G = self.copy()
        G.remove_nodes_from([n for n in G if n not in other])
        return G

    def draw(self, prog="dot", attribute="fillcolor",  hold=False, **kargs):
        """Draw the network using matplotlib. Not exactly what we want but could be useful.

        :param str prog: one of the graphviz program (default dot)
        :param bool hold: hold previous plot (default is False)
        :param str attribute: attribute to use to color the nodes (default is "fillcolor").
        :param node_size: default 1200
        :param width: default 2

        Uses the fillcolor attribute of the nodes
        Uses the link attribute of the edges

        .. seealso:: :meth:`plot` that is dedicated to this kind of plot using graphviz


        """
        self.logging.warning("Not for production. Use plot() instead")
        pos = nx.drawing.graphviz_layout(self, prog=prog)

        node_size = kargs.get('node_size', 1200)
        kargs['node_size'] = node_size

        width = kargs.get('width', 2)
        kargs['width'] = width

        # node attributes
        nodes = sorted(self.nodes())
        node_colors = [self.node[x][attribute] if attribute in self.node[x].keys()
                else "gray" for x in nodes]

        # edge attributes
        edges = self.edges(data=True)
        colors = {'-':'red', '+':'black'}
        edge_colors = [colors[x[2]['link']] for x in edges]

        nx.draw(self, prog=prog, hold=hold, nodelist=nodes,
            edge_color=edge_colors, node_color=node_colors,
            pos=pos, **kargs)

    def _check_dot_prog(self, prog):
        DevTools().check_param_in_list(prog, ["twopi", "gvcolor", "wc", "ccomps", "tred",
            "sccmap", "fdp", "circo", "neato", "acyclic", "nop", "gvpr", "dot",
            "sfdp"])

    def _get_cmap(self, cmap=None):
        if cmap == "heat":
            cmap = self._colormap.get_cmap_heat_r()
        elif cmap == "green":
            cmap = self._colormap.get_cmap_red_green()
        else:
            try:
                cmap = colormap.cmap_builder(cmap)
            except:
                pass # already a valid cmap ?
        return cmap

    def plot(self, prog="dot", viewer="pylab", hold=False,
        show=True, filename=None, node_attribute=None, edge_attribute=None,
        cmap='heat', colorbar=False, remove_dot=True,
        normalise_cmap=True, edge_attribute_labels=True,
        rank_method='inout'
        ):
        """plotting graph using dot program (graphviz) and networkx

        By default, a temporary file is created to hold the image created by
        graphviz, which is them shown using pylab. You can choose not to see the
        image (show=False) and to save it in a local file instead (set the
        filename). The output format is PNG. You can play with
        networkx.write_dot to save the dot and create the SVG yourself.

        :param str prog: the graphviz layout algorithm (default is dot)
        :param viewer: pylab
        :param bool show: show the plot (True by default)
        :param bool remove_dot: if True, remove the temporary dot file.
        :param edge_attribute_labels: is True, if the label are available, show them.
            otherwise, if edge_attribute is provided, set lael as the edge_attribute

        :param rank_method: If none, let graphviz do the job. Issue is that (i)
            input stimuli may not be aligned and output neither. The rank_method set
            to **cno** constraints the stimuli and measured species that are sinks
            all others are free. The **same** constraint all nodes with same rank
            to be in the same subgraph.

        Additional attributes on the graph itself can be set up by populating the
        graph attribute with a dictionary called "graph"::

            c.graph['graph'] = {"splines":True, "size":(20,20), "dpi":200}

        Useful other options are::

            c.edge["tnfa"]["tnfr"]["penwidth"] = 3
            c.edge["tnfa"]["tnfr"]["label"] = " 5"

        If you use edge_attribute and show_edge_labels, label are replaced
        by the content of edge_attribute. If you still want differnt labels,
        you must stet show_label_edge to False and set the label attribute
        manually

        ::

            c = cnograph.CNOGraph()
            c.add_reaction("A=C")
            c.add_reaction("B=C")
            c.edge['A']['C']['measure'] = 0.5
            c.edge['B']['C']['measure'] = 0.1
            c.expand_and_gates()
            c.edge['A']['C']['label'] = "0.5 seconds"
            # compare this that shows only one user-defined label
            c.plot()
            # with that show all labels
            c.plot(edge_attribute="whatever", edge_attribute_labels=False)

        See the graphviz homepage documentation for more options.

        ::

            c.plot(filename='test.svg', viewer='yout_favorite_viewer',
                remove_dot=False, rank_method='cno')


        .. note:: edge attribute in CNOGraph (Directed Graph) are not coded
            in the same way in CNOGraphMultiEdges (Multi Directed Graph).
            So, this function does not work for MultiGraph

        .. todo:: use same colorbar as in midas. rigtht now; the vmax is not correct.
        .. todo:: precision on edge_attribute to 2 digits.

        if filename provided with extension different from png, pylab must be able to
        read the image. If not, you should set viewer to something else.



        c.plot() # uses default layout, with PNG format viewed with pylab
        c.plot(rank_method='same') # uses ranks of each node with PNG format viewed with pylab
        c.plot(rank_method='cno', filename='test.svg', viewer='browse') # uses browse (defautl viewer)
        c.plot(rank_method='cno', filename='test.svg', viewer='browse')

        """
        # graph is a DiGraph attribute
        # that is sometimes replaced by {} inside networkx so we need to overwrite it here
        # each time we want to plot the graph.
        if len(self) == 0:
            self.logging.error("empty graph, nothing to plot")
            return

        if rank_method:
            assert rank_method in ['inout', 'all', 'cno']
        self._check_dot_prog(prog)

        # Get the default/user attributes for the graph/nodes/edges for graphviz
        # FIXME is this still required ?
        #self.graph = self.graph_options.copy()

        # Set the colormap
        cmap = self._get_cmap(cmap)

        # update the node attributes if required with default color
        # required is we manipualte the _signals, _inbitors and so on
        self.set_default_node_attributes()
        # or ues the requried node attribute.
        M = 1
        if node_attribute is not None:
            # TODO check that it exists
            # cmap = matplotlib.cm.get_cmap(cmap)
            sm = matplotlib.cm.ScalarMappable(
                norm = matplotlib.colors.Normalize(vmin=0, vmax=1), cmap=cmap)

            # node[0] is the name, node[1] is the data
            data = [node[1][node_attribute] for node in self.nodes(data=True)
                    if node_attribute in node[1].keys()]
            if normalise_cmap is True:
                M = max(data)

            # color could be encoded as values between 0 and 1
            # or hexa. Default to 1. If not all are provided,
            # no errors raised.
            for node in self.nodes():
                # default
                self.node[node]['fillcolor'] = "#FFFFFF"
                try:
                    value = self.node[node][node_attribute]/float(M)
                    rgb = sm.to_rgba(value)
                    colorHex = matplotlib.colors.rgb2hex(rgb)
                    self.node[node]['fillcolor'] = colorHex
                except:
                    try:
                        color = self.node[node][node_attribute]
                        self.node[node]['fillcolor'] = colormap.Color(color).hex
                    except:
                        pass

        # update the edge attribute
        if edge_attribute:
            M = self._set_edge_attribute_color(edge_attribute, cmap)

            # set edge with edge_attribute set to 0 as invisible edges
            reactions = [self.edge2reaction(edge) for edge in self.edges(data=True)
                    if edge[2][edge_attribute] > 0]
            self.set_edge_visibility_from_reactions(reactions)


        # create temp files
        # FIXME we create png here ?? do we use outfile ?
        infile  = tempfile.NamedTemporaryFile(suffix=".dot", delete=False)
        if filename is None:
            outfile = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
            filename = outfile.name

        # Some nodes may belong to 2 colors. Creating subgraph is one way to go
        # around. Having graphviz 2.30 we could use striped node.
        if node_attribute is None:
            for node in self.nodes():
                if node in self.signals and node in self.inhibitors:
                    self.node[node]['style'] = "diagonals,filled"
                    self.node[node]['color'] = "red"
                    self.node[node]['fillcolor'] = "lightblue"
                if node in self.stimuli and node in self.inhibitors:
                    self.node[node]['style'] = "diagonals,filled"
                    self.node[node]['color'] = "red"
                    self.node[node]['fillcolor'] = "#9ACD32"

        # to not change the current graph, let us copy it
        # FIXME we use 'this' variable  and use it for edges
        # but not for the nodes... why ?
        this = self.copy()

        if edge_attribute_labels and edge_attribute is not None:
            self._set_edge_attribute_label(this, edge_attribute)

        count = 0
        ret = -1
        while count < 10:
            # this is a hack for graphviz 2.30
            H = self._get_ranked_agraph(rank_method)
            H.write(infile.name)
            frmt = os.path.splitext(filename)[1][1:]
            try:
                if filename.endswith('svg') and 'dpi' in H.graph_attr.keys():
                    del H.graph_attr['dpi']
                if rank_method != 'cno':
                    H.draw(path=filename, prog=prog, format=frmt)
                else:
                    H.to_dot(infile.name, frmt=frmt)
                    args=[" -T"+frmt, infile.name]
                    args=' '.join(args)
                    data = H._run_prog(prog, args)
                    fh = open(filename, "w")
                    fh.write(data)
                    fh.close()

                count = 1000
                ret = 0
            except Exception as err:
                print(err.message)
                self.logging.warning("%s program failed. Trying again" % prog)
                count += 1

        if ret !=0:
            if rank_method is not None:
                self.logging.warning("%s program failed to create image" % prog)
            H = self._get_ranked_agraph('cno')
            #H = nx.to_agraph(this)
            if filename.endswith('svg') and 'dpi' in H.graph_attr.keys():
                del H.graph_attr['dpi']
            frmt = os.path.splitext(filename)[1][1:]
            #H.draw(path=filename, prog=prog, format=frmt)
            H.to_dot(infile.name, frmt=frmt)
            args=[" -T"+frmt, infile.name]
            args=' '.join(args)
            #print(os.path.exists(infile.name))
            data = H._run_prog(prog, args)
            fh = open(filename, "w")
            fh.write(data)

        # Here is the visualisation only iif image was created
        if viewer=="pylab" and show is True and frmt!='svg':
            if hold == False:
                if hold is False:
                    pylab.clf()
                f = pylab.gcf()
                f.set_facecolor("white")
                a = f.add_axes([0.05,0.04,0.8,0.9])
                a.imshow(pylab.imread(filename))
                a.axis('off')
            else:
                a = pylab.imshow(pylab.imread(filename))
                pylab.axis('off')

            if colorbar:
                cbar = pylab.linspace(0, 1, 100)
                e = f.add_axes([.86, .05, .05, .9])
                e.pcolor(pylab.array([cbar, cbar]).transpose(),
                        cmap=cmap, vmin=0, vmax=1);
                e.yaxis.tick_right()
                e.set_xticks([],[])

                e.set_yticks([0,20,40,60,80,100])
                if normalise_cmap is True:
                    # TODO precision for small number ??
                    from easydev.tools import precision
                    e.set_yticklabels([precision(x,3) for x in pylab.linspace(0, M, 6)])
                else:
                    e.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])

        elif show is True:
            try:
                if viewer=='pylab':
                    viewer = 'browse'
                subprocess.call("%s %s &" % (viewer, filename), shell=True)
            except:
                # for MAC users ?
                if viewer!='pylab':
                    subprocess.call("open -a Preview %s &" % (viewer, filename), shell=True)

        if remove_dot == True:
            infile.delete = True
            infile.close()
        else:
            print("Creating dot file named 'model.dot'")
            shutil.move(infile.name, "model.dot")
            infile.close()

        if filename == False:
            outfile.delete = True
            outfile.close()

        if node_attribute is not None:
            self.set_default_node_attributes()

        #if show == True:
        #    try:
        #        from biokit.dev.mpl_focus import ZoomPan
        #        ax = pylab.gca()
        #        zp = ZoomPan()
        #        _ = zp.zoom_factory(ax, base_scale = 1.2)
        #        _ = zp.pan_factory(ax)
        #    except:
        #        pass

    def _repr_png_(self):
        """Returns an Image for display in an IPython console"""
        fh = TempFile(suffix='.png')
        self.plot(show=False, filename=fh.name, 
                rank_method=self.graph_options['ipython']['rank_method'])
        return fh.name

    def _repr_svg_(self):
        """Returns an Image for display in an IPython console"""
        from IPython.display import SVG
        fh = TempFile(suffix='.svg')
        self.plot(show=False, filename=fh.name)
        return fh.name

    @property
    def png(self):
        from IPython.display import Image
        data = Image(self._repr_png_(), embed=True, 
                width=self.graph_options['ipython']['width'])
        return data

    @property
    def svg(self):
        from IPython.display import SVG
        data = SVG(self._repr_svg_() )
        return data

    @property
    def _ipython_display_(self):
        from IPython.display import Image
        data = Image(self._repr_png_(), embed=True, width=200)
        return data

    def _plot_rmse_fromR(self, filename, F=.3, scale=2, col=None):
        """
        filename should be a CSV with first column being node names and second
        column the values.

        This function populates the node attribute "mse", which should be the
        average RMSE over all experiment for each species.

        :param col: if col is None, average all columns

        c.plot()
        """
        import pandas as pd
        df = pd.read_csv(filename, index_col=0,  header=None)

        if col==None:
            df = df.mean(axis=1)
        elif col in df.columns:
            df = df[col]
        else:
            print("Invalid column provided. Use one of {}".format(df.columns))
        for this in self.nodes():
            if this in self.signals:
                mse = df.ix[this] #.values[0]
                self.node[this]['mse'] =  (1-(mse/F)**scale)
                self.node[this]['label'] =  this+"\n"+str(int(mse*1000)/1000.)
            else:
                self.node[this]['mse'] = 1
        cm = colormap.Colormap()
        self.plot(node_attribute="mse", cmap=cm.get_cmap_heat())

    def _set_edge_attribute_label(self, this, edge_attribute):
        for e in this.edges():
            this.edge[e[0]][e[1]]['label'] = this.edge[e[0]][e[1]][edge_attribute]

    def _set_edge_attribute_color(self, edge_attribute, cmap):
        import matplotlib
        cmap = matplotlib.cm.get_cmap(cmap)
        sm = matplotlib.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(vmin=0, vmax=1), cmap=cmap)

        # normalise if needed only
        values = [self.edge[edge[0]][edge[1]][edge_attribute] for edge in self.edges()]
        if max(values)>1:
            M = max([self.edge[edge[0]][edge[1]][edge_attribute] for edge in self.edges()])
        else:
            M = 1.

        for edge in self.edges():
            value = self.edge[edge[0]][edge[1]][edge_attribute] / M
            rgb = sm.to_rgba(value)
            colorHex = matplotlib.colors.rgb2hex(rgb)
            self.edge[edge[0]][edge[1]]['color'] = colorHex
        return M

    def _get_ranked_agraph(self, rank_method=None):
        """and gates should have intermediate ranks"""

        if rank_method == 'cno':
            H = AGraphCNO(self)
        else:
            H = nx.to_agraph(self)

        for k, v in self.graph_options['graph'].items():
            H.graph_attr[k] = v
        for k, v in self.graph_options['edge'].items():
            H.edge_attr[k] = v
        for k, v in self.graph_options['node'].items():
            H.node_attr[k] = v

        if rank_method is None:
            return H
        if self.midas is None:
            return H
        if rank_method == 'cno':
            return H
        # order the graph for ranks
        allranks = self.get_same_rank() # this has been checkd on MMB to and seems correct
        ranks  = {}
        M = max(allranks.keys())
        for k, v in allranks.items():
            if rank_method in ['inout', 'all']:
                #H.strict = True
                ranks[k] = sorted([x for x in v if '=' not in x],
                        key=lambda x: x.lower())
            else:
                ranks[k] = [x for x in v if '=' not in x]

        selfloops = [x[1] for x in self.selfloop_edges()]

        # Note: if name is set to "cluster"+name, black box is added
        # Group by cluster with same ranks

        # Note2 that  subgraph should be added before the invisiable edges.
        # The inverse resets the attribute e.g., style="invis"

        for rank in ranks.keys():
            name = str(rank) # may be numbers
            if rank == 0:
                # label will be used if name == 'cluster_source'
                species = [x for x in ranks[rank] if x not in selfloops]
                aa = H.add_subgraph(species,  rank='source', name='source',
                        label='stimuli')
                if len(species)>=2:
                    for i, node1 in enumerate(species[0:-1]):
                        node2 = species[i+1]
                        if self._graph_type != 'multigraph':
                            aa.add_edge(node1, node2, style="invis")
            elif rank == M:
                species = sorted([x for x in ranks[rank] if x not in selfloops])
                aa = H.add_subgraph(species, name="sink", rank='sink')
                if len(species)>=2:
                    for i, node1 in enumerate(species[0:-1]):
                        node2 = species[i+1]
                        if self._graph_type != 'multigraph':
                            aa.add_edge(node1, node2, style="invis")
            else:
                if rank_method in ["all", 'cnor']:
                    # something funky here with self loop
                    # self loop should be ignored here
                    species = sorted([x for x in ranks[rank] if x not in selfloops])
                    H.add_subgraph(species, name=name, rank='all')

        return H

    def _get_nonc(self):
        if self._nonc is None:
            nonc = self.findnonc()
            self._nonc = nonc
        return self._nonc
    nonc = property(fget=_get_nonc,
        doc="Returns list of Non observable and non controlable nodes (Read-only).")

    def _get_reactions(self):
        # todo: could use edge2reactions
        # tood: sort reactions inside SIF
        sif = self.to_sif()
        return sorted(sif.reactions)
    reactions = property(_get_reactions, doc="return the reactions (edges)")

    def _get_namesSpecies(self):
        nodes = self.nodes()
        nodes = [x for x in nodes if "+" not in x and "=" not in x]
        return sorted(nodes)
    species = property(fget=_get_namesSpecies,
        doc="Return sorted list of species (ignoring and gates) Read-only attribute.")

    #@do_profile()
    def swap_edges(self, nswap=1, inplace=True, self_loop=False):
        """Swap two edges in the graph while keeping the node degrees fixed.

        A double-edge swap two randomly chosen edges u-v and x-y
        and creates the new edges u-x and v-y::

            u  v                u  v
            |  |     becomes    |  |
            x  y                y  x

        If either the edge u-x or v-y already exist no swap is performed
        and another attempt is made to find a suitable edge pair.

        :param int nswap: number of swaps to perform (Defaults to 1)
        :return: nothing

        .. warning:: the graph is modified in place.

        .. warning:: and gates are removed at the beginning since they do not make sense
            with the new topolgy. They can easily be added back 

        a proposal swap is ignored in 3 cases:
        #. if the summation of in_degree is changed
        #. if the summation of out_degree is changed
        #. if resulting graph is disconnected

        .. note:: what about self loop ? if proposed, there are ignored except 
            if required to be kept

        """
        self._changed = True
        self.remove_and_gates()
        Ninh = [x[2]["link"] for x in self.edges(data=True)].count('-')
        I = sum(self.in_degree().values())
        O = sum(self.out_degree().values())

        # find 2 nodes that have at least one successor
        count = 0
        trials = 0
        status = {'and':0, 'selfloop':0, 'indegree':0, 'outdegree':0, 'unconnected':0}
        if inplace is False:
            local = self.copy()

        while count < nswap and trials<nswap*5:
            trials += 1
            if inplace is False:
                edges = local.edges()
            else:
                edges = self.edges()
            np.random.shuffle(edges)
            e1, e2 = edges[0:2]

            # ignore self loop:
            if (e1[0] == e2[1] or e2[0] == e1[1]) and self_loop is False:
                status['selfloop']+=1
                continue

            if inplace is False:
                d1 = local.edge[e1[0]][e1[1]].copy()
                d2 = local.edge[e2[0]][e2[1]].copy()
                G = nx.DiGraph(local)
            else:
                d1 = self.edge[e1[0]][e1[1]].copy()
                d2 = self.edge[e2[0]][e2[1]].copy()
                G = nx.DiGraph(self)

            #self.copy()
            G.add_edge(e1[0], e2[1], None, **d1)
            G.add_edge(e2[0], e1[1], None, **d2)
            G.remove_edge(e1[0], e1[1])
            G.remove_edge(e2[0], e2[1])
            #print(count, len(G))

            if sum(G.in_degree().values()) != I:
                # the link already exists
                status['indegree'] +=1
                continue

            # This is slow, let us use something different
            #if nx.is_connected(G.to_undirected()) == False:
            #    status['unconnected'] += 1
            #    continue
            if 0 in G.degree().values():
                status['unconnected'] += 1
                continue

            if sum(G.out_degree().values()) != O:
                status['outdegree'] +=1
                continue

            # seems okay , so let us swap the edges now.
            if inplace is True:
                self.add_edge(e1[0], e2[1], None, **d1)
                self.add_edge(e2[0], e1[1], None, **d2)
                self.remove_edge(e1[0], e1[1])
                self.remove_edge(e2[0], e2[1])
            else:
                local.add_edge(e1[0], e2[1], None, **d1)
                local.add_edge(e2[0], e1[1], None, **d2)
                local.remove_edge(e1[0], e1[1])
                local.remove_edge(e2[0], e2[1])

            # number of inhibitory link must remain identical
            #Ninh2 = [x[2]["link"] for x in self.edges(data=True)].count('-')
            #assert Ninh2 == Ninh
            # seems to be true 
            #assert nx.is_connected(self.to_undirected()) == True
            count +=1
        status['count'] = count
        status['trials'] = trials

        self._status_swap_edges = status

        if inplace is False:
            return local, status
        else:
            self._changed = True

    def adjacency_matrix(self, nodelist=None, weight=None):
        """Return adjacency matrix.

        :param list nodelist: The rows and columns are ordered according to the nodes in nodelist.
            If nodelist is None, then the ordering is produced by :meth:`nodes` method.

        :param str weight: (default=None) The edge data key used to provide each value in the matrix.
            If None, then each edge has weight 1. Otherwise, you can set it to
            "weight"

        :returns: numpy matrix Adjacency matrix representation of CNOGraph.

        .. note:: alias to :meth:`networkx.adjacency_matrix`

        .. seealso:: :meth:`adjacency_iter` and :meth:`adjacency_list`

        """
        return nx.adjacency_matrix(self, nodelist=nodelist).astype(int)

    def remove_edge(self, u, v):
        """Remove the edge between u and v.

        :param str u: node u
        :param str u: node v

        Calls :meth:`clean_orphan_ands` afterwards
        """
        self._changed = True
        super(CNOGraph, self).remove_edge(u,v)
        self.clean_orphan_ands()

    @modifier
    def remove_node(self, n):
        """Remove a node n

        :param str node: the node to be removed

        Edges linked to **n** are also removed. **AND** gate may now be
        irrelevant (only one input or no input) and are therefore removed as well.

        .. seealso:: :meth:`clean_orphan_ands`

        """
        self._changed = True
        super(CNOGraph, self).remove_node(n)
        if "^" not in str(n):
            self.clean_orphan_ands()

    def add_nodes_from(self, nbunch, attr_dict=None, **attr):
        self._changed = True
        for node in nbunch:
            if isinstance(node, str):
                self.add_node(node, attr_dict=attr_dict, **attr)
            elif isinstance(node, tuple):
                name = node[0]
                attr.update(node[1])
                self.add_node(name, attr_dict=attr_dict, **attr)

    def add_node(self, node, attr_dict=None, **attr):
        """Add a node

        :param str node: a node to add
        :param dict attr_dict: dictionary, if None,
            replace by default shape/color. recognised keys are those
            from graphviz.
            Default keys set are color, fillcolor, shape, style, penwidth
        :param attr: additional keyword provided will overwrite keys found in
            **attr_dict** parameter

        ::

            c = CNOGraph()
            c.add_node("A", data=[1,2,3,])
            c.add_node("B", shape='circle')

        .. warning:: **attr** replaces any key found in attr_dict. See :meth:`add_edge` for details.

        .. todo:: currently nodes that contains a ^ sign are interpreted as AND gate and will appear
           as small circle. One way to go around is to use the label attribute.
           you first add the node with a differnt name and populate the label with
           the correct nale (the one that contain the ^ sign); When calling the plot
           function, they should all appear as expected.

        """
        self._changed = True
        if attr_dict is None:
            attr_dict = self.get_node_attributes(node)
        super(CNOGraph, self).add_node(node, attr_dict=attr_dict, **attr)

    def preprocessing(self, expansion=True, compression=True, cutnonc=True,
                      maxInputsPerGate=2):
        """Performs the 3 preprocessing steps (cutnonc, expansion, compression)

        :param bool expansion: calls :meth:`expand_and_gates` method
        :param bool compression: calls :meth:`compress` method
        :param bool cutnon: calls :meth:`cutnonc` method
        :param int maxInputPerGates: parameter for the expansion step

        .. plot::
            :width: 80%
            :include-source:

            from cno import CNOGraph, cnodata
            c = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
            c.preprocessing()
            c.plot()

        """
        self._changed = True
        if cutnonc:
            self.cutnonc()
        if compression:
            self.compress()
        if expansion:
            self.expand_and_gates(maxInputsPerGate=maxInputsPerGate)

    def cutnonc(self):
        """Finds non-observable and non-controllable nodes and removes them.

        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph, cnodata
            c = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
            c.cutnonc()
            c.plot()

        """
        self._changed = True
        nonc = self.nonc[:]
        for node in nonc:
            self.collapse_node(node)

    def compress(self, recursive=True, iteration=1, max_iteration=5):
        """Finds compressable nodes and removes them from the graph

        A compressable node is a node that is not part of the special nodes
        (stimuli/inhibitors/readouts mentionned in the MIDAS file). Nodes
        that have multiple inputs and multiple outputs are not compressable
        either.


        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph, cnodata
            c = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
            c.cutnonc()
            c.compress()
            c.plot()



        .. seealso:: :meth:`compressable_nodes`, :meth:`is_compressable`
        """
        self._changed = True
        assert max_iteration >=1
        #Patched to proceed on sorted lists and provide always the same results.
        #for node in self.compressable_nodes:
        for node in sorted(self.compressable_nodes):
            if self.is_compressable(node) == True:
                # update the graph G 
                self._compressed.append(node)
                self.collapse_node(node)

        # proceed on a sorted list and provide 
        # always the same results in theory.
        #
        for node in sorted(self.nodes()):
            if self.degree(node) == 0 and node not in self.stimuli and \
                node not in self.signals and node not in self.inhibitors:
                """
                FIXME It consists in finding the list of nodes that belong
                to the network but are in no edges (if a node A is
                removed during compression, it appears in no edges,
                but you can still find it in CNOGraph.nodes()). For
                this set of nodes that do not belong to the network,
                I add fake edges like A + mock1 and so on.

                Doing like this I am able to save the compressed (or
                the repaired) network, while using the same MIDAS file.

                I was afraid to touch the code because I am not sure
                whether this is the intended behaviour."""
                self.logging.info("Found an orphan, which has been removed (%s)" % node)
                self.remove_node(node)

        if len(self.compressable_nodes) > 0:
            self.logging.warning("There are still compressable nodes. Call again")
            if recursive and iteration<max_iteration:
                self.compress(iteration=iteration+1)

        # finally, let us rename the AND gates if needed
        self._update_and_gate_names()

    def _get_and_attribute(self, reaction):
        # attributes on the N edges of an AND gate must be identical
        preds = self.predecessors(reaction)
        attr = self.edge[preds[0]][reaction]
        return attr

    def _update_and_gate_names(self):
        """Update and gates
        
        names are updated. Individual names are sorted alphabetically
        """
        ands_before = self._find_and_nodes()
        #
        for node in ands_before:
            preds = self.signed_predecessors(node)
            succs = self.successors(node)
            if len(succs) == 1:
                succ = succs[0]
                # build new AND gate and sort it alphabetically
                newname = Reaction("=".join(["^".join(preds), succ]))
                newname.sort()
                # sort the current node alphabetically as well.
                r1 = Reaction(node)
                r1.sort()
                
                if newname != r1.name:
                    attr = self._get_and_attribute(node)
                    attr2 = self.attributes.attributes['and']
                    self.add_reaction(newname.name)
                    self.node[newname.name].update(**attr2)
                    self.node[newname.name].update(**attr)
                    self.remove_node(node)
            else:
                CNOError("an AND gate cannot have several successors. This should not happen")

    def signed_predecessors(self, node):
        preds = self.predecessors(node)
        for i, pred in enumerate(preds):
            if self.edge[pred][node]['link'] == '-':
                preds[i] = "!" + pred
        return preds

    def _get_collapse_edge(self, inputAttrs, outputAttrs):
        attrs = inputAttrs.copy()
        if inputAttrs['link'] == outputAttrs['link']:
            attrs['link'] = "+"
            attrs['color'] = "black"
            attrs['arrowhead'] = "normal"
        else:
            attrs['link'] = "-"
            attrs['color'] = "red"
            attrs['arrowhead'] = "tee"
        return attrs

    def collapse_node(self, node):
        """Collapses a node (removes a node but connects input nodes to output nodes)

        This is different from :meth:`remove_node`, which removes a node and its edges thus
        creating non-connected graph. :meth:`collapse_node`, instead remove the node but merge the input/output
        edges IF possible. If there are multiple inputs AND multiple outputs the
        node is not removed.

        :param str node: a node to collapse.

        * Nodes are collapsed if there is at least one input or output.
        * Node are not removed if there is several inputs and several ouputs.
        * if the input edge is -, and the next is + or viceversa then the final edge if -
        * if the input edge is - and output is - then final edge is +

        """
        # Removing a node that has several entries and several outputs is not
        # implemented because there is no such situation in CNO.
        self._changed = True
        successors = self.successors(node)
        predecessors = self.predecessors(node)

        # todo: may be simplified ?
        if len(successors) == 1 and len(predecessors)==1:
            self.logging.debug("Compressing %s 1,1 mode" % node)
            # compressed node
            attr1 = self.edge[node][successors[0]]
            #
            attr2 = self.edge[predecessors[0]][node]

            attrs = self._get_collapse_edge(attr1, attr2)
            if predecessors[0] != successors[0]:
                self.add_edge(predecessors[0], successors[0], None, **attrs)
        elif len(successors) == 1:

            for predecessor in predecessors:
                attr = self.edge[predecessor][node]
                if predecessor != successors[0]:
                    attr2 = self.edge[node][successors[0]]
                    attrs = self._get_collapse_edge(attr, attr2)
                    self.add_edge(predecessor, successors[0], None, **attrs)
        elif len(predecessors) == 1:
            for successor in successors:
                attr = self.edge[node][successor]
                if predecessors[0] != successor:
                    attr2 = self.edge[predecessors[0]][node]
                    attrs = self._get_collapse_edge(attr, attr2)
                    self.add_edge(predecessors[0], successor, None,  **attrs)
        else:
            if len(successors) > 1 and len(predecessors) > 1:
                self.logging.debug(node, successors, predecessors)
                self.logging.warning("N succ >1 and N pred >1. Node not removed. use remove_node() if you really want to remove it")
                return
                #raise ValueError("invalid node to remove several in/out")
            else:
                self.logging.debug("%s %s %s" % (node, successors, predecessors))
                self.logging.warning("unknown case (no output or input ?). Node %s removed"% node)
                #raise ValueError("invalid node to remove several in/out")
        self.remove_node(node)
        self._changed = True

    def get_node_attributes(self, node):
        """Returns attributes of a node using the MIDAS attribute

        Given a node, this function identifies the type of the input
        node and returns a dictionary with the relevant attributes found
        in :attr:`node_attributes.attributes`.

        For instance, if a midas file exists and if **node** belongs to the stimuli,
        then the dicitonary returned contains the color green.

        :param str node:
        :returns: dictionary of attributes.

        """
        # default
        attr = self.attributes['others'].copy()
        # otherwisen
        if self.midas:
            if node in self.stimuli:
                attr = self.attributes['stimuli'].copy()
            elif node in self.signals:
                attr = self.attributes['signals'].copy()
            elif node in self.inhibitors:
                attr = self.attributes['inhibitors'].copy()
            elif node in self.nonc: # or node in self.findnonc():
                attr = self.attributes['nonc'].copy()
            elif '^' in node:
                attr = self.attributes['and'].copy()
            elif node in self.compressable_nodes:
                attr = self.attributes['compressed'].copy()
        else:
            if '^' in str(node):
                attr = self.attributes['and'].copy()

        if node in self._stimuli:
            attr = self.attributes['stimuli'].copy()
        if node in self._signals:
            attr = self.attributes['signals'].copy()
        if node in self._inhibitors:
            attr = self.attributes['inhibitors'].copy()
        return attr

    def predecessors_as_reactions(self, node):
        """Build up reactions for a given node from predecessors only"""
        predecessors = self.predecessors(node)
        reactions = []
        for pred in predecessors:
            if self.isand(pred):
                reactions.append(pred)
            elif self[pred][node]['link'] == '+':
                reactions.append(pred + "=" +node)
            else:
                reactions.append("!" + pred + "=" +node)
        return reactions

    def set_default_node_attributes(self):
        """Set all node attributes to default attributes

        .. seealso:: :meth:`get_node_attributes`
        """
        if self._set_default_node_attribute == 'cno':
            for node in self.nodes():
                attrs = self.get_node_attributes(node)
                for k, v in attrs.items():
                    self.node[node][k] = v

    def set_node_attributes(self, attr_dict ):
        """Experimental

        dictionaries where keys are nodes and keys are
        dictionaries of attributes.  Does not need to
        be

        {'Akt': {'color':red', 'strength':0.5}}
        Does not need to have a value for each node
        could set default values may be ?
        """

        for node in attr_dict.keys():
            if node in self.nodes():
                for k,v in attr_dict[node].items():
                    self.node[node][k] = v

    def _remove_edge_attribute(self, name):
        for edge in self.edges(data=True):
            try:
                del self.edge[edge[0]][edge[1]][name]
            except:
                pass


    def set_edge_attribute(self, name, values):
        # FIXME: not used can be deleted
        # Used in plot_optimised_model in base.py
        # values is a dictionary of reaction as keys
        for reaction in values.keys():
            if reaction in values.keys():
                edges = self.reac2edges(reaction)
                for edge in edges:
                    e0, e1 = edge[0], edge[1]
                    self.edge[e0][e1][name] = values[reaction]

    def reac2edges(self, reaction):
        """Here no informatiom about links is returned"""
        #if reaction not in self.reactions:
        #    raise CNOError("Unknown reaction {0}".format(reaction))
        reac = Reaction(reaction)
        dic = reac.get_signed_lhs_species()

        if "^" in reaction:
            inputs = reac.lhs_species
            out = reac.rhs

            def get_links(this):
                if this in dic['+']:
                    return '+'
                else:
                    return '-'
            return [(this, reaction, get_links(this)) for this in inputs] + [(reaction,out)]
        else: #simple case
            if dic['+']:
                link = '+'
            else:
                link ='-'

            e0 = reac.lhs
            e1 = reac.rhs
            return [(e0.replace("!",""), e1.replace("!",""), link)]

    def set_edge_visibility_from_reactions(self, reactions):
        # First, reset the style, which may have been set already
        self._remove_edge_attribute('style')

        # AND gate should be handled differently

        for edge in self.edges(data=True):
            if self.and_symbol in edge[1]:
                reaction = edge[0] + "=" + edge[1]
                if edge[2]['link'] == "-":
                    reaction = "!" + reaction
            elif self.and_symbol in edge[0]:
                reaction = edge[0] + "=" + edge[1]
                if edge[2]['link'] == "-":
                    reaction = "!" + reaction
            else:
                reaction = self.edge2reaction(edge)

            if self._graph_type != 'multigraph':
                if reaction not in reactions:
                    self.edge[edge[0]][edge[1]]['style'] = 'invis'

    def get_same_rank(self):
        """Return ranks of the nodes.


        """
        # no need to run this function, which could be long if already computed
        if self._ranks is not None:
            return self._ranks.copy()

        import time
        t1 = time.time()
        # some aliases
        try:
            stimuli = self.stimuli
            if len(stimuli) == 0:
                stimuli = self._get_inputs()
        except:
            stimuli = self._get_inputs()

        try:
            signals = self.signals
            if len(signals) == 0:
                signals = self._get_outputs()
        except:
            signals = self._get_outputs()

        func_path = nx.algorithms.floyd_warshall(self)

        maxrank = int(self.get_max_rank())
        # start populating the ranks starting with the obvious one: stimuli and
        # signals
        ranks = collections.defaultdict(list)
        ranks[0] = stimuli
        for node in sorted(self.nodes(),
                key=lambda x: str(x).lower()):
            # skip and gate
            if self.isand(node):
                continue
            # skip end signals for now
            if node in signals and len(self.successors(node))==0:
                continue
            elif node not in stimuli:
                distances = [func_path[s][node] for s in stimuli]
                distances = [x for x in distances if x != pylab.inf]
                if len(distances) != 0:
                    M = int(np.nanmax([abs(x) for x in distances if x != pylab.inf]))
                    ranks[M].append(node)
                else:
                    self.logging.debug('warning, rank %s is empyt'% node)

        # now the end signal
        for node in sorted(self.nodes(), key=lambda x: str(x).lower()):
            if node in signals and len(self.successors(node))==0:
                ranks[maxrank+1].append(node)

        t2 = time.time()
        if t2-t1 > 1:
            print("get_same_rank %s seconds" % str(t2-t1))
        self._ranks = ranks.copy()
        return ranks

    def _get_inputs(self):
        return [x for x in self.nodes() if len(self.predecessors(x))==0]

    def _get_outputs(self):
        self.logging.warning('WARNING. no signals found, tryring to build a list (node with no successors)')
        return  [x for x in self.nodes() if len(self.successors(x))==0]

    def get_max_rank(self):
        """Get the maximum rank from the inputs using floyd warshall algorithm

        If a MIDAS file is provided, the inputs correspond to the stimuli.
        Otherwise, (or if there is no stimuli in the MIDAS file),
        use the nodes that have no predecessors as inputs (ie, rank=0).

        """
        stimuli = self.stimuli
        if len(stimuli) == 0:
            stimuli = self._get_inputs()

        func_path = nx.algorithms.floyd_warshall(self)
        # compute the longest path from Stimuli by using the floyd warshall
        # algorithm. inputs/stimuli has rank 0.
        ranks = [[x for x in func_path[stimulus].values() if x !=pylab.inf]
            for stimulus in stimuli]

        allranks = []
        for r in ranks:
            allranks = allranks + r #concatenate all ranks includeing empty []
        maxrank = np.nanmax(allranks)
        return maxrank

    def _add_and_gates(self, node, maxInputsPerGate=2):
        """See expand_and_gates docstring"""
        preds = self.predecessors(node)
        preds = [pred for pred in preds if self.isand(pred) is False]
        assert maxInputsPerGate>=2 and maxInputsPerGate<=5, "maxInputsPerGate must be >2 and less than 5"
        #todo: order predecessirs according to nameSpecies order
        self.logging.debug( "node %s, pred=%s " % (node,  preds))

        if len(preds) == 1:
            self.logging.debug("Nothing to do with %s (only 1 predecessor)" % node)
            return

        for inputsPerGates in range(2, maxInputsPerGate+1):
            self.logging.debug(inputsPerGates)
            if inputsPerGates>len(preds):
                continue
            for combi in itertools.combinations(preds, inputsPerGates):
                self.logging.debug("adding %s input gate from %s to %s" % (inputsPerGates, combi, node))
                shape = "circle"
                if len(combi) == 3:
                    shape = "triangle"
                elif len(combi) >= 4:
                    shape = "square"
                andNode = self._nodes2reac(list(combi), node)
                self.logging.debug("add_and_gates: %s " % andNode)
                self.add_node(andNode, shape=shape)
                for combinode in combi:
                    attr = self.edge[combinode][node].copy()
                    attr['link'] = self.edge[combinode][node]['link']
                    #attr['weight'] = pylab.nan # edge from and to specy always + and nan weight
                    self.add_edge(combinode, andNode, None, **attr)
                attr['link'] = "+" # edge from and to specy always +
                #attr['weight'] = pylab.nan # edge from and to specy always + and nan weight
                attr['color'] = 'black' # and output is always black
                self.add_edge(andNode, node, None, **attr)
        self._changed = True

    def _nodes2reac(self, inputsNodes, output):
        inputs = []
        for node in inputsNodes:
            if self.edge[node][output]['link']=="-":
                inputs.append("!"+str(node))
            else:
                inputs.append(node)

        reac = self.and_symbol.join([str(x) for x in inputs])
        reac += "=" + str(output)
        reac = Reaction(reac)
        reac.sort()
        return reac.name

    def edge2reaction(self, edge):
        e1 = edge[0]
        e2 = edge[1]
        link = edge[2]['link']
        if link == '+':
            return "=".join([e1, e2])
        else:
            return "=".join(["!"+e1, e2])

    def isand(self, node):
        if self.and_symbol in str(node):
            return True
        else:
            return False

    def _find_nodes_with_multiple_inputs(self):
        """return a list of nodes that have multiple predecessors"""
        nodes = []
        for node in self.nodes():
            if len(self.predecessors(node)) > 1 and self.isand(node) is False:
                nodes.append(node)
            else:
                if len(self.predecessors(node)) > 1 and self.isand(node):
                    self.logging.debug("ignore " + str(node))
        return nodes

    def _find_and_nodes(self):
        andNodes = [node for node in self.nodes() if self.isand(node)]
        return andNodes

    def expand_or_gates(self):
        """Expand OR gates given AND gates

        If a graph contains AND gates (without its OR gates), you can add back
        the OR gates automatically using this function.

        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph
            from pylab import subplot, title

            c1 = CNOGraph()
            c1.add_edge("A", "C", link="-")
            c1.add_edge("B", "C", link="+")
            c1.expand_and_gates()
            subplot(1,3,1)
            title("OR and AND gates")
            c1.plot(hold=True)

            c1.remove_edge("A", "C")
            c1.remove_edge("B", "C")
            subplot(1,3,2)
            c1.plot(hold=True)
            title("AND gates only")

            c1.expand_or_gates()
            subplot(1,3,3)
            c1.plot(hold=True)
            title("after call to \\n expand_or_gates function")

        .. seealso:: :meth:`~cno.io.cnograph.CNOGraph.expand_and_gates`

        """
        for this in self._find_and_nodes():
            p = self.predecessors(this)
            s = self.successors(this)
            assert len(s) == 1
            for node in p:
                link = self.edge[node][this]['link']
                self.add_edge(node, s[0], link=link)
        self._changed = True

    def expand_and_gates(self, maxInputsPerGate=2):
        """Expands the network to incorporate AND gates

        :param int maxInputsPerGate: restrict maximum number of inputs used to
            create AND gates (default is 2)

        The CNOGraph instance can be used to model a boolean network.  If a node
        has several inputs,  then the combinaison of the inputs behaves like an OR gate
        that is we can take the minimum over the inputs.

        In order to include AND behaviour, we introduce a special node called
        AND gate. This function adds AND gates whenever a node has several
        inputs. The AND gates can later on be used in a boolean formalism.

        In order to recognise AND gates, we name them according to the following
        rule. If a node A has two inputs B and C, then the AND gate is named::

            B^C=A

        and 3 edges are added: B to the AND gates, C to the AND gates and the AND gate to A.

        If an edge is a "-" link then, an **!** character is introduced.

        In this expansion process, AND gates themselves are ignored.

        If there are more than 2 inputs, all combinaison of inputs may be
        considered but the default parameter **maxInputsPerGate** is set to 2.
        For instance, with 3 inputs A,B,C you may have the
        following combinaison: A^B, A^C, B^C. The link A^B^C will be added only
        if **maxInputsPerGate** is set to 3.


        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph
            from pylab import subplot, title

            c = CNOGraph()
            c.add_edge("A", "C", link="+")
            c.add_edge("B", "C", link="+")
            subplot(1,2,1)
            title("Original network")
            c.plot(hold=True)

            c.expand_and_gates()
            subplot(1,2,2)
            c.plot(hold=True)
            title("Expanded network")

        .. seealso:: :meth:`remove_and_gates`, :meth:`clean_orphan_ands`,
            :meth:`expand_or_gates`.

        .. note:: this method adds all AND gates in one go. If you want to add a specific AND gate,
            you have to do it manually. You can use the :meth:`add_reaction` for that purpose.


        .. note:: propagate data from edge on the AND gates.
        """
        self._changed = True
        nodes2expand = self._find_nodes_with_multiple_inputs()
        for node in nodes2expand:
            self._add_and_gates(node, maxInputsPerGate)

    def add_cycle(self, nodes, **attr):
        """Add a cycle

        :param list nodes: a list of nodes. A cycle will be constructed from
           the nodes (in order) and added to the graph.
        :param dict attr: must provide the "link" keyword. Valid values are "+", "-"
            the links of every edge in the cycle will be identical.

        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph
            c = CNOGraph()
            c.add_edge("A", "C", link="+")
            c.add_edge("B", "C", link="+")
            c.add_cycle(["B", "C", "D"], link="-")
            c.plot()

        .. warning:: added cycle overwrite previous edges

        """
        self._changed = True
        if "link" not in attr.keys():
            raise KeyError("link keyword must be provided")

        attr = self.get_default_edge_attributes(**attr)
        super(CNOGraph, self).add_cycle(nodes, **attr)

    def add_path(self):
        """networkx method not to be used"""
        raise NotImplementedError

    def add_star(self):
        """networkx method not to be used"""
        raise NotImplementedError

    def remove_edges_from(self, edges):
        for edge in edges:
            self.remove_edge(edge)

    def add_weighted_edges_from(self):
        """networkx method not to be used"""
        raise NotImplementedError


    def remove_nodes_from(self, nbunch):
        """Removes a bunch of nodes

        .. warning:: need to be tests with and gates."""
        self._changed = True
        for node in nbunch:
            self.remove_node(node)

    def _get_compressable_nodes(self):
        if self._changed is True:
            compressables = [x for x in self.nodes() if self.is_compressable(x)]
            self._compressable_nodes = compressables[:]
            self._changed = False
        else:
            compressables = self._compressable_nodes[:]
        return [x for x in compressables if self.isand(x) is False]
    compressable_nodes = property(fget=_get_compressable_nodes,
        doc="Returns list of compressable nodes (Read-only).")

    def _ambiguous_multiedge(self, node):
        """Test whether the removal of a node may lead to multi edges ambiguity

        e.g., A-->B and A--| B
        """
        edges = [(x[0], x[1], x[2]['link']) for x in self.edges(data=True)]
        for A in self.predecessors(node):
            for B in self.successors(node):
                link = self.edge[node][B]['link']
                if link == "+":
                    link = "-"
                elif link == "-":
                    link = "+"
                if (A, B, link) in edges:
                    return True
        return False

    def is_compressable(self, node):
        """Returns True if the node can be compressed, False otherwise

        :param str node: a valid node name
        :return: boolean value

        Here are the rules for compression. The main idea is that a node can be removed if the
        boolean logic is preserved (i.e. truth table on remaining nodes is preserved).

        A node is compressable if it is not part of the stimuli, inhibitors, or
        signals specified in the MIDAS file.

        If a node has several outputs and inputs, it cannot be compressed.

        If a node has one input or one output, it may be compressed.
        However, we must check the following possible ambiguity that could be raised by the removal
        of the node: once removed, the output of the node may have multiple input edges with different
        types of inputs edges that has a truth table different from the original truth table.
        In such case, the node cannot be compressed.

        Finally, a node cannot be compressed if one input is also an output (e.g., cycle).


        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph
            from pylab import subplot,show, title
            c = cnograph.CNOGraph()
            c.add_edge("a", "c", link="-")
            c.add_edge("b", "c", link="+")
            c.add_edge("c", "d", link="+")
            c.add_edge("b", "d", link="-")
            c.add_edge("d", "e", link="-")
            c.add_edge("e", "g", link="+")
            c.add_edge("g", "h", link="+")
            c.add_edge("h", "g", link="+")

            # multiple inputs/outputs are not removed
            c.add_edge("s1", "n1", link="+")
            c.add_edge("s2", "n1", link="+")
            c.add_edge("n1", "o1", link="+")
            c.add_edge("n1", "o2", link="+")

            c._stimuli = ["a", "b", "s1", "s2"]
            c._signals = ["d", "g", "o1", "o2"]

            subplot(1,2,1)
            c.plot(hold=True)
            title("Initial graph")

            c.compress()
            subplot(1,2,2)
            c.plot(hold=True)
            title("compressed graph")

            show()


        """
        if node not in self.nodes():
            msg = "node %s is not in the graph" % node
            raise ValueError(msg)

        # there is always a MIDAS file for now but we may remove it later on
        specialNodes = self.inhibitors + self.signals  + self.stimuli
        notcompressable = [x for x in self.nodes() if x in specialNodes]
        notcompressable += [x for x in self.nodes() if self.isand(x)]
        if node in notcompressable:
            return False

        succs = self.successors(node)
        preds = self.predecessors(node)

        # if one input is an input, no compression
        if len(set(succs).intersection(set(preds))) > 0:
            self.logging.debug('skipped node (retroaction ? %s) ' % node)
            #print("%s case cycle input is in output" % node)
            return False
        # if input and ouputs are provided, they should contain not AND on
        # both the input and output
        # see https://github.com/cellnopt/cellnopt/issues/39
        a = sum([1 for x in succs if self.isand(x)])
        b = sum([1 for x in preds if self.isand(x)])
        if a>=1 and b>=1:
            return False

        # if no output or no input, can be compressed.
        # somehow this is redundant with NONC algorithm, which is not required
        # anymore.
        if len(self.successors(node)) == 0 or len(self.predecessors(node))==0:
            self.logging.debug('skipped node %s (input/output case ' % node)
            return True

        # if multiple input AND multiple output, no ambiguity: nothing to compress
        elif len(self.successors(node)) >1 and len(self.predecessors(node))>1:
            self.logging.debug('skipped node %s >1,>1 case ' % node)
            return False
        # if one input and several outputs OR several input and one ouptut, we
        # can compress. However, if the compressable node once removed is an
        # edge that already exists, then we will have multiedges, which is not
        # handled in a DIgraph. So, we cannot comrpess that particular case.

        # if one input only and one output only, no ambiguity: we can compress
        elif len(self.predecessors(node)) == 1 and len(self.successors(node))==1:
            ambiguous = self._ambiguous_multiedge(node)
            if ambiguous == True:
                self.logging.debug('%s could be compressed but ambiguity with existing edge so not removed ' % node)
                return False
            else:
                self.logging.debug('Add node %s =1,=1 to be removed ' % node)
                return True

        # one output but several input may be ambiguous. If output is inhibitor
        # and input contains an mix of inhibitor/activation then we can not
        # compress
        elif len(preds) > 1 and len(succs)==1:
            input_links = [self[p][node]['link'] for p in preds]
            output_links = self[node][succs[0]]['link']
            if (output_links == "-") and len(set(input_links))>=2:
                self.logging.debug('skipped node %s output=inhibitor and ambiguous inputs' % node)
                return False
            else:
                ambiguous = self._ambiguous_multiedge(node)
                if ambiguous == True:
                    self.logging.debug('%s could be compressed but ambiguity with existing edge so not removed ' % node)
                    return False
                else:
                    self.logging.debug('Add node %s >1,=1 case to be removed' % node)
                    return True

        # one input and several ouptut, no ambiguity: we can compress
        elif len(preds) == 1 and len(succs)>1:
            ambiguous = self._ambiguous_multiedge(node)
            if ambiguous == True:
                self.logging.debug('%s could be compressed but ambiguity with existing edge so not removed ' % node)
                return False
            else:
                self.logging.debug('Add node %s =1,>1 case to be removed' % node)
                return True
        else:
            self.logging.debug('do not remove node %s' % node)
            return False

    def relabel_nodes(self, mapping, inplace=True):
        """Function to rename a node, while keeping all its attributes.

        :param dict mapping: a dictionary mapping old names (keys) to new names
            (values )
        :return: new :class:`CNOGraph` object


        if we take this example::

            c = CNOGraph();
            c.add_reaction("a=b");
            c.add_reaction("a=c");
            c.add_reaction("b=d");
            c.add_reaction("c=d");
            c.expand_and_gates()

        Here, an AND gate has been created. c.nodes() tells us that its name is
        "b^c=d". If we rename the node b to blong, the AND gate name is
        unchanged if we use the nx.relabel_nodes function. Visually, it is
        correct but internally, the "b^c=d" has no more meaning since the node
        "b" does not exist anymore. This may lead to further issues if we for
        instance split the node c::

            c = nx.relabel_nodes(c, {"b": "blong"})
            c.split_node("c", ["c1", "c2"])

        This function calls relabel_node taking care of the AND nodes as well.

        .. warning:: inplace modifications

        .. todo:: midas must also be modified

        """
        # keep names of the AND gates for book keeping
        # ands_before = self._find_and_nodes()
        # mapping_ands = {}
        for this in self._find_and_nodes():
            reac = Reaction(this)
            reac.rename_species(mapping)
            # add to mapping
            if reac.name != this:
                mapping.update({this:reac.name})

        self.mapping = mapping

        # ideally, nx.relabels_nodes(self, mapping, copy=False) could
        # do the work but somehow loses the edge data dictionary.,,
        # so, we will do it ourself:
        edge_data = [x for x in self.edges(data=True) if x[0] in mapping.keys()
            or x[1] in mapping.keys()]

        self.edge_data = edge_data
        for edge in edge_data:
            # do not remove the edge first otherwise AND may be lost
            # since clean_orphans may be called.

            # now add new edges. First,
            if "=" in edge[0]:  # an AND gate
                r = Reaction(edge[0])
                r.rename_species(mapping)
                e0 = r.name
            elif edge[0] in mapping.keys():
                e0 = mapping[edge[0]]
            else:
                # normal node
                e0 = edge[0]

            if "=" in edge[1]:  # an AND gate
                r = Reaction(edge[1])
                r.rename_species(mapping)
                e1 = r.name
            elif edge[1] in mapping.keys():
                e1 = mapping[edge[1]]
            else:
                e1 = edge[1]
            self.add_edge(e0, e1, **edge[2])

        # Now, handle the MIDAS file if any by renaming all species updating
        # the stimuli/inhibitors/signals if needed. Need to keep track of user
        # defined names
        self._stimuli = [mapping[x] if x in mapping.keys() else x for x in self._stimuli]
        self._signals = [mapping[x] if x in mapping.keys() else x for x in self._signals]
        self._inhibitors = [mapping[x] if x in mapping.keys() else x for x in self._inhibitors]
        self._compressed = [mapping[x] if x in mapping.keys() else x for x in self._compressed]

        # copy the node data
        for k, v in mapping.items():
            self.node[v] = self.node[k]

        if self.midas is not None:
            for k, v in mapping.items():
                if "^" in k:
                    del mapping[k]
            try: self.midas.rename_species(mapping)
            except: pass
            try: self.midas.rename_inhibitors(mapping)
            except: pass
            try: self.midas.rename_stimuli(mapping)
            except: pass

        # finally remove the old nodes
        self.remove_nodes_from(mapping.keys())
        self._changed = True

    def centrality_eigenvector(self, max_iter=1000, tol=0.1):
        res = nx.eigenvector_centrality(self,max_iter,tol=tol)
        nx.set_node_attributes(self, 'eigenvector', res)
        nx.set_node_attributes(self, 'centrality_eigenvector', res)
        import operator
        degcent_sorted = sorted(res.items(), key=operator.itemgetter(1), reverse=True)
        #for k,v in degcent_sorted:
        #    self.logging.info("Highest degree centrality %s %s", v,k)
        return res

    def centrality_degree(self):
        """Compute the degree centrality for nodes.

        The degree centrality for a node v is the fraction of nodes it
        is connected to.

        :return: list of nodes with their degree centrality. It is also added to the list
            of attributes with the name "degree_centr"

        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph, cnodata
            c = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
            c.centrality_degree()
            c.plot(node_attribute="centrality_degree")


        """
        res = nx.degree_centrality(self)
        nx.set_node_attributes(self, 'degree', res)
        nx.set_node_attributes(self, 'centrality_degree', res)
        import operator
        degcent_sorted = sorted(res.items(), key=operator.itemgetter(1), reverse=True)
        #for k,v in degcent_sorted:
        #    self.logging.info("Highest degree centrality %s %s", v,k)
        return res

    def centrality_closeness(self, **kargs):
        """Compute closeness centrality for nodes.

        Closeness centrality at a node is 1/average distance to all other nodes.

        :param v: node, optional  Return only the value for node v
        :param str distance: string key, optional (default=None)  Use specified edge key as edge distance.   If True, use 'weight' as the edge key.
        :param bool normalized:  optional   If True (default) normalize by the graph size.

        :return: Dictionary of nodes with closeness centrality as the value.

        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph, cnodata
            c = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
            c.centrality_closeness()
            c.plot(node_attribute="centrality_closeness")

        """
        res = nx.centrality.closeness_centrality(self, **kargs)

        nx.set_node_attributes(self, 'closeness', res)
        nx.set_node_attributes(self, 'centrality_closeness', res)
        import operator
        degcent_sorted = sorted(res.items(), key=operator.itemgetter(1), reverse=True)
        #for k,v in degcent_sorted:
        #    self.logging.info("Highest closeness centrality %s %s", v,k)
        return res

    def centrality_betweeness(self, k=None, normalized=True,
        weight=None, endpoints=False, seed=None):
        r"""Compute the shortest-path betweeness centrality for nodes.

        Betweenness centrality of a node `v` is the sum of the
        fraction of all-pairs shortest paths that pass through `v`:

        .. math::

            c_B(v) =\sum_{s,t \in V} \frac{\sigma(s, t|v)}{\sigma(s, t)}

        where :math:`V` is the set of nodes, :math:`\sigma(s, t)` is the number of
        shortest :math:`(s, t)`-paths,  and :math:`\sigma(s, t|v)` is the number of those
        paths  passing through some  node :math:`v` other than :math:`s, t`.
        If :math:`s = t`, :math:`\sigma(s, t) = 1`, and if :math:`v \in {s, t}`,
        :math:`\sigma(s, t|v) = 0` .

        :param int k: (default=None)
          If k is not None use k node samples to estimate betweeness.
          The value of k <= n where n is the number of nodes in the graph.
          Higher values give better approximation.
        :param bool normalized:   If True the betweeness values are normalized by :math:`2/((n-1)(n-2))`
          for graphs, and :math:`1/((n-1)(n-2))` for directed graphs where :math:`n`
          is the number of nodes in G.
        :param str weight: None or string, optional
          If None, all edge weights are considered equal.

        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph, cnodata
            c = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
            c.centrality_betweeness()
            c.plot(node_attribute="centrality_betweeness")

        .. seealso:: networkx.centrality.centrality_betweeness

        """
        res = nx.centrality.betweenness_centrality(self,k=k,normalized=normalized,
            weight=weight, endpoints=endpoints, seed=seed)
        nx.set_node_attributes(self, 'centrality_betweeness', res)
        nx.set_node_attributes(self, 'betweeness', res)
        import operator
        degcent_sorted = sorted(res.items(), key=operator.itemgetter(1), reverse=True)
        #for k,v in degcent_sorted:
        #    self.logging.info("Highest betweeness centrality %s %s", v,k)
        return res

    def to_gexf(self, filename):
        """Export into GEXF format

        :param str filename:

        This is the networkx implementation and requires the version 1.7
        This format is quite rich and can be used in external software such as Gephi.

        .. warning:: color and labels are lost. information is stored as
            attributes.and should be as properties somehow.
            Examples:  c.node['mkk7']['viz'] =  {'color': {'a': 0.6, 'r': 239, 'b': 66,'g': 173}}

        """
        from networkx.readwrite import write_gexf
        write_gexf(self, filename)

    def to_directed(self):
        """networkx method not to be used"""
        raise NotImplementedError

    def to_sif(self, filename=None):
        """Export CNOGraph into a SIF file.

        Takes into account and gates. If a species called  "A^B=C" is found, it is an AND
        gate that is encoded in a CSV file as::

            A 1 and1
            B 1 and1
            and1 1 C

        :param str filename:

        """
        sif = SIF()

        for edge in self.edges(data=True):
            n1 = edge[0]
            n2 = edge[1]
            link = edge[2]['link']
            reaction = ""
            if self.isand(n1) is False and self.isand(n2) is False:
                if link == "-":
                    reaction += "!"
                reaction += n1 + "=" + n2
                sif.add_reaction(reaction)
        for edge in self._find_and_nodes():
            sif.add_reaction(edge)

        if filename:
            sif.save(filename)
        else:
            return sif

    def to_json(self, filename=None):
        """Export the graph into a JSON format

        :param str filename:

        .. seealso:: :meth:`loadjson`
        """
        from networkx.readwrite import json_graph
        data = json_graph.node_link_data(self)
        if filename is not None:
            json.dump(data, open(filename, "w"))
        else:
            return data

    def read_sbmlqual(self, filename):
        sif = SIF()
        sif.read_sbmlqual(filename)
        self.clear()
        for reac in sif.reactions:
            self.add_reaction(reac)

    def to_sbmlqual(self, filename=None):
        """Export the topology into SBMLqual and save in a file

        :param str filename: if not provided, returns the SBML as a string.
        :return: nothing if filename is not provided

        .. seealso:: :meth:`cno.io.sbmlqual`
        """
        self._changed = True
        s = SIF()
        for reac in self.reactions:
            s.add_reaction(reac)
        return s.to_sbmlqual(filename=filename)

    def read_json(self, filename):
        """Load a network in JSON format as exported from :meth:`to_json`

        :param str filename:

        .. seealso:: :meth:`to_json`
        """
        self._changed = True
        from networkx.readwrite import json_graph

        graph = json_graph.node_link_graph(json.loads(open(filename).read()))
        self.clear()

        for node in graph.nodes():
            self.add_node(node)
        for edge in graph.edges(data=True):
            self.add_edge(edge[0], edge[1], link=edge[2]['link'])

        return graph

    def lookfor(self, specyName):
        """Prints information about a species

        If not found, try to find species by ignoring cases.

        """
        #try to find the specy first:
        if specyName not in self.nodes():
            print("did not find the requested specy",)
            proposals = [node for node in self.nodes() if specyName.lower() in node.lower()]
            if len(proposals):
                print("try one of ")
                for p in proposals:
                    print(proposals)
        else:
            print("predecessors")
            print(self.predecessors(specyName))
            print("successors")
            print(self.successors(specyName))

    def get_stats(self):
        stats = {}
        flow = nx.flow_hierarchy(self)
        stats['flow'] = flow
        stats['mean_degree'] = sum(self.degree().values())/float(len(self.nodes()))
        return stats

    def summary(self):
        """Plot information about the graph"""
        stats = self.get_stats()
        print("Flow hierarchy = %s (fraction of edges not participating in cycles)" % stats['flow'])
        print("Average degree = " + str(sum(self.degree().values())/float(len(self.nodes()))))

    def merge_nodes(self, nodes, node):
        """Merge several nodes into a single one

        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph
            from pylab import subplot
            c = CNOGraph()
            c.add_edge("AKT2", "B", link="+")
            c.add_edge("AKT1", "B", link="+")
            c.add_edge("A", "AKT2", link="+")
            c.add_edge("A", "AKT1", link="+")
            c.add_edge("C", "AKT1", link="+")
            subplot(1,2,1)
            c.plot(hold=True)
            c.merge_nodes(["AKT1", "AKT2"], "AKT")
            subplot(1,2,2)
            c.plot(hold=True)


        """
        self._changed = True
        assert len(nodes)>1, "nodes must be a list of species"
        for n in nodes:
            if n not in self.nodes():
                raise ValueError("%s not found in the graph !" % n)

        self.add_node(node)
        for n in nodes:
            for pred in self.predecessors(n):
                attrs = self.edge[pred][n]
                if self.isand(pred):
                    pred = self._rename_node_in_reaction(pred, n, node)
                    self.add_reaction(pred)
                else:
                    self.add_edge(pred, node, **attrs)
            for succ in self.successors(n):
                attrs = self.edge[n][succ]
                if self.isand(succ):
                    succ = self._rename_node_in_reaction(succ, n, node)
                    self.add_reaction(succ)
                else:
                    self.add_edge(node, succ, **attrs)
        for n in nodes:
            self.remove_node(n)

        # assume that all nodes are in signals if first one is in signal.
        # FIXME: not robust
        if nodes[0] in self._signals:
            for this in nodes:
                self._signals.remove(this)
            self._signals.append(node)
        if nodes[0] in self._stimuli:
            for this in nodes:
                self._stimuli.remove(this)
            self._stimuli.append(node)
        if nodes[0] in self._inhibitors:
            for this in nodes:
                self._inhibitors.remove(this)
            self._inhibitors.append(node)

        self._changed = True

    def split_node(self, node, nodes):
        """

        .. plot::
            :include-source:
            :width: 50%

            from cno import CNOGraph
            from pylab import subplot
            c = CNOGraph()
            c.add_reaction("!A=C")
            c.add_reaction("B=C")
            c.add_reaction("!b1=B")
            c.add_reaction("b2=B")
            c.expand_and_gates()

            subplot(1,2,1)
            c.plot(hold=True)

            c.split_node("B", ["B1", "B2", "B3"])
            subplot(1,2,2)
            c.plot(hold=True)

        """
        for n in nodes:
            for pred in self.predecessors(node):
                attrs = self.edge[pred][node]
                if self.isand(pred):
                    pred = self._rename_node_in_reaction(pred, node, n)
                    self.add_reaction(pred)
                else:
                    self.add_edge(pred, n, **attrs)

            for succ in self.successors(node):
                attrs = self.edge[node][succ]
                # special case of the AND gates
                if self.isand(succ):
                    succ = self._rename_node_in_reaction(succ, node, n)
                    self.add_reaction(succ)
                else: # normal case
                    self.add_edge(n, succ, **attrs)
        self.remove_node(node)
        # remove AND gates as well:
        for this in self._find_and_nodes():
            gate = Reaction(this)
            if node in gate.lhs_species:
                self.remove_node(this)

        if node in self._signals:
            self._signals.extend(nodes)
            self._signals.remove(node)
        if node in self._stimuli:
            self._stimuli.extend(nodes)
            self._stimuli.remove(node)
        if node in self._inhibitors:
            self._inhibitors.extend(nodes)
            self._inhibitors.remove(node)
        self._changed = True

    def _rename_node_in_reaction(self, reaction, old, new):
        """This function rename a species within a reaction."""
        reac = Reaction(reaction)
        reac.rename_species({old:new})
        return reac.name

    def findnonc(self):
        """Finds the Non-Observable and Non-Controllable nodes

        #. Non observable nodes are those that do not have a path to any measured
           species in the PKN
        #. Non controllable nodes are those that do not receive any information
           from a species that is perturbed in the data.

        Such nodes can be removed without affecting the readouts.


        :param G: a CNOGraph object
        :param stimuli: list of stimuli
        :param stimuli: list of signals

        :return: a list of names found in G that are NONC nodes

        .. doctest::

            >>> from cno import CNOGraph, cnodata
            >>> model = cnodata('PKN-ToyMMB.sif')
            >>> data = cnodata('MD-ToyMMB.csv')
            >>> c = CNOGraph(model, data)
            >>> namesNONC = c.nonc()

        :Details: Using a floyd Warshall algorithm to compute path between nodes in a
          directed graph, this class
          identifies the nodes that are not connected to any signals (Non Observable)
          and/or any stimuli (Non Controllable) excluding the signals and stimuli,
          which are kept whatever is the outcome of the FW algorithm.
        """
        # some aliases
        assert (self.stimuli!=None and self.signals!=None)

        dist = nx.algorithms.floyd_warshall(self)

        namesNONC = []
        for node in self.nodes():
            # search for paths from the species to the signals
            spe2sig = [(node, dist[node][s]) for s in self.signals if dist[node][s]!=np.inf]
            # and from the nstimuli to the species
            sti2spe = [(node, dist[s][node]) for s in self.stimuli if dist[s][node]!=np.inf]

            if len(spe2sig)==0 or len(sti2spe)==0:
                if node not in self.signals and node not in self.stimuli:
                    namesNONC.append(node)

        namesNONC  = list(set(namesNONC)) # required ?
        return namesNONC

    def remove_self_loops(self):
        self._changed = True
        for e in self.edges():
            if e[0] == e[1]:
                self.remove_edge(e[0], e[1])



class AGraphCNO(gv.AGraph):

    def __init__(self, N):
        """N should be a cnograph so that we can get the ranks"""
        self.cnograph = N

        directed = N.is_directed()
        strict = N.number_of_selfloops()==0 and not N.is_multigraph()

        super(AGraphCNO, self).__init__(name=N.name, strict=strict, directed=directed)


        # default graph attributes
        self.graph_attr.update(N.graph.get('graph',{}))
        self.node_attr.update(N.graph.get('node',{}))
        self.edge_attr.update(N.graph.get('edge',{}))

        # add nodes
        for n,nodedata in N.nodes(data=True):
            self.add_node(n,**nodedata)

        # loop over edges

        if N.is_multigraph():
            for u,v,key,edgedata in N.edges_iter(data=True,keys=True):
                str_edgedata=dict((k,str(v)) for k,v in edgedata.items())
                self.add_edge(u,v,key=str(key),**str_edgedata)
        else:
            for u,v,edgedata in N.edges_iter(data=True):
                str_edgedata=dict((k,str(v)) for k,v in edgedata.items())
                self.add_edge(u,v,**str_edgedata)
        self._changed = True

    def to_dot(self, filename=None, frmt='png'):
        allranks = self.cnograph.get_same_rank()

        txt = "digraph G{\n"
        for k, v in self.cnograph.graph_options['graph'].items():
            if 'svg' in frmt and k == 'dpi':
                continue
            txt += """    %s="%s";\n""" % (k,v)

        for rank in sorted(allranks.keys()):
            if rank == 0:
                rankname = "source"
            elif rank == max(allranks.keys()):
                rankname = "sink"
            else:
                rankname = 'same'
            if len(allranks[rank]):
                names = ";".join(sorted(allranks[rank]))
                txt += "    {rank=%s;%s;}\n" % (rankname, names)

        for node in self.cnograph.nodes(data=True):
            attrs = node[1].copy()
            if "^" in node[0]:
                attrs['name'] = '"%s"' % node[0]
                attrs['label'] = ''
                attrs['penwidth'] = '' + 'shape="circle" width=.1 height=.1 fixedsize=True'
            else:
                attrs['name'] = node[0]
                attrs['label'] = node[0]
                try:
                    attrs['penwidth'] = """penwidth="%s" """ % attrs['penwidth']
                except:
                    print("FIXME")
                    attrs['penwidth'] = ''
                    pass
            txt += """%(name)s [color="%(color)s", fillcolor="%(fillcolor)s" %(penwidth)s  shape="%(shape)s" style="%(style)s" label="%(label)s" fontname=Helvetica fontsize=22.0 ];\n""" % attrs

        for edge in self.cnograph.edges(data=True):
            attrs = edge[2]
            attrs['in'] = '"%s"' % edge[0]
            attrs['out'] = '"%s"' % edge[1]
            attrs['label']= ''
            txt += """ %(in)s -> %(out)s [ color="%(color)s" label="%(label)s" weight="1.000000" penwidth="2" arrowhead="%(arrowhead)s" style="solid"];\n""" % (attrs)
        txt += "}"

        fh = open(filename, 'w')
        fh.write(txt)
        fh.close()

