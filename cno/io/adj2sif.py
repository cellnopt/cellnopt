# -*- python -*-
#
#  This file is part of cellnopt software
#
#  Copyright (c) 2012-2014 - EMBL-EBI
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>, 
#      <cokelaer at gmail dot com>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: www.cellnopt.org
#
##############################################################################
""":Topic: **adjacency matrix**"""
import networkx as nx
import numpy


class ADJ2SIF(object):
    """Reads an adjacency matrix (and names) from CSV files

    .. warning:: API likely to change to use pandas to simplify the API.

    The instance can then be exported to :class:`~cno.io.sif.SIF` or used
    as input for the :class:`cno.io.cnograph.CNOGraph` structure.

    ::

        >>> from cno.io import ADJ2SIF
        >>> from cno import getdata
        >>> f1 = getdata("test_adjacency_matrix.csv")
        >>> f2 = getdata("test_adjacency_names.csv")
        >>> s = ADJ2SIF(f1, f2)
        >>> sif = s.to_sif()
        >>> c = CNOGraph(s.G)
    
        Where the adjacency matrix looks like::

            0,1,0
            1,0,0
            0,0,1

        and names is a 1-column file::

            A
            B
            C

        The exported SIF file would look like::

            A 1 B
            A 1 C

    .. warning:: The adjacency matrix contains only ones (no -1)  so future
        version may need to add that information using incidence matrix for
        instance


    """
    def __init__(self, filenamePKN=None, filenameNames=None, delimiter=","):
        """.. rubric:: Constructor

        :param str filenamePKN: adjacency matrix made of 0's and 1's. 
        :param str filenameNames: names of the columns/rows of the adjacency matrix
        :param str delimiter: commas by default

        ::

            0,1,0
            1,0,0
            0,0,1

        names::

            A
            B
            C

        The 2 files above correspond to this SIF file::

            A 1 B
            A 1 C

        """
        self.filename = filenamePKN
        self.filenameNames = filenameNames
        self.delimiter = delimiter

        self._G = None
        self._names = None

        if self.filename:
            self.load_adjacency()

        if filenameNames:
            self.load_names()

    def _get_names(self):
        return self._names
    names = property(_get_names, 
        doc="Names of the nodes read from the the provided filename. Could be empty")

    def _get_G(self):
        return self._G
    G = property(_get_G, doc="The graph created from the input data")

    def load_adjacency(self, filename=None):
        """Reads an adjacency matrix filename
        
        if no filename is provided, tries to load from the attribute
        :attr:`filename`.
        
        
        """
        if filename:
            self.filename = filename
        self._G = nx.Graph(numpy.loadtxt(self.filename, delimiter=self.delimiter))

    def load_names(self, filename=None):
        """Reads the columns/rows names"""
        if filename:
            self.filenameNames = filename

        fh = open(self.filenameNames, "r")
        data = fh.read()
        fh.close()
        self._names = data.split()

 
    def to_sif(self, filename=None):
        """Exports input data files into a SIF instance and save it

        :param str filename: set this parameter if you want to save the SIF into
            a file
        :return: a SIF instance

        """
        from cno.io.sif import SIF
        s = SIF()
        for edge in self.G.edges_iter():
            reac = self.names[edge[0]] + "=" + self.names[edge[1]]
            s.add_reaction(reac)
        if filename: 
            s.save(filename)
        return s


