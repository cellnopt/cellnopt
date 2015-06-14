# -*- python -*-
#
#  This file is part of cellnopt.core software
#
#  Copyright (c) 2011-2014 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: www.cellnopt.org
#
##############################################################################
import csv
from .sif import SIF


class EDA(object):
    """Reads networks in EDA format

    EDA format is similar to SIF but provides a weight on each edge.

    It looks like::

        A (1) B = .5
        B (1) C =  1
        A (1) C = .1

    .. note:: the parentheses and spaces.
    .. note:: no header expected.
    """
    def __init__(self, filename, threshold=0, verbose=False):
        """

        :param str filename:
        :param float threshold: should be between 0 and 1 but not compulsary
        :param bool verbose:


        """
        self.filename = filename
        self.verbose = verbose
        self.threshold = threshold
        try:
            self._load()
        except Exception as err:
            print(err.message)
            print("Could not read the EDA file")
            raise Exception

    def _load(self):
        if self.verbose:
            print("Loading EDA scores from %s" % self.filename),
        data = open(self.filename)
        data_iter = csv.reader(data, delimiter=" ")
        data_iter = list(data_iter)

        if len(data_iter) and "edgescore" in [x.lower() for x in data_iter[0]]:
            del data_iter[0]

        # skip first row
        #header = data_iter.next()
        data = []
        for i, datum in enumerate(data_iter):
            if len(datum) == 5:
                data.append(datum)
            elif len(datum) == 0 or len(datum) == 1:
                print("Found empty line in EDA file %s" %  self.filename)
            else:
                self.error("- Line %s in %s is ill formed (expected 5 columns): %s" % (i,self.filename, datum))
                break
        self.nodes1 = [x[0] for x in data]
        self.nodes2 = [x[2] for x in data]
        self.edges = [x[1] for x in data]
        self.scores = [float(x[4]) for x in data]

    def export2sif(self, threshold=None):
        """Exports EDA data into SIF file

        :param float threshold: since EDA format provides a weight on each edge,
            it can be used as a threshold to consider the relation or not.
            By default, the :attr:`threshold` is set to 0, which means all edges
            should be exported in the output SIF format (assuming weights are positive).
            You ca n either set the :attr:`threshold` attribute to a different value
            or provide this **threshold** parameter to override the default threshold.

        ::

            >>> from cno.io.eda import EDA
            >>> from cno import testing
            >>> e = EDA(testing.get("test_simple.eda"))
            >>> s1 = e.export2sif() # default threshold 0
            >>> len(s1)
            3
            >>> s1 = e.export2sif(0.6) # one edge with weight=0.5 is ignored
            >>> len(s1)
            2

        """
        if threshold == None:
            threshold = self.threshold

        s = SIF()
        for n1, edge, score, n2 in zip(self.nodes1, self.edges, self.scores, self.nodes2):
            if score > threshold:
                if edge == '(1)':
                    s.add_reaction("%s=%s" % (n1, n2))
                elif edge == '(-1)':
                    s.add_reaction("!%s=%s" % (n1, n2))

        return s

    def plot(self):
        """

        c.add_edge('E', 'F', label='0.1', width=2, penwidth=2, link="+")
        c.plot(edge_attribute='penwidth', cmap='gray_r')

        """
        from cno import CNOGraph
        c = CNOGraph()
        for n1, edge, score, n2 in zip(self.nodes1, self.edges, self.scores, self.nodes2):
            if edge == '(1)':
                c.add_edge(n1, n2, link='+', width=score, penwidth=score)
            elif edge == '(-1)':
                c.add_edge(n1, n2, link='-', width=score, penwidth=score)
        c.plot(edge_attribute='penwidth', cmap='gray_r')


