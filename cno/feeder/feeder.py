# -*- python -*-
#
#  This file is part of cellnopt software
#
#  Copyright (c) 2014 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: http://github.com/cellnopt
#
##############################################################################
"""Wrapping of CNORfeeder package (see bioconductor)"""

from cno.core import CNORBase


class Feeder(CNORBase):
    """Find missing links automatically based on CNORfeeder R package

    .. plot::
        :include-source:
        :width: 80%

        from cno import Feeder, cnodata
        feeder = Feeder()
        feeder.run(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"),
                   k=4, verbose=False)
        feeder.newlinks
        feeder.plot()
    """
    def __init__(self, verboseR=True):
        CNORBase.__init__(self, verboseR)

        self.Rscript = """
            library(CNORfeeder)
            pknmodel = readSIF("%(pkn)s")
            cnolist = CNOlist("%(data)s")

            model = preprocessing(cnolist, pknmodel,
                                  compression=%(compression)s,
                                  expansion=%(expansion)s,
                                  maxInputsPerGate=%(maxInputsPerGate)s)
            # find new links
            BTable = makeBTables(CNOlist=cnolist, k=%(k)s,
                                measErr=c(%(err1)s,%(err2)s))
            alllinks = linksRanking(CNOlist=cnolist,
                                    measErr=c(%(err1)s,%(err2)s))
            modelIntegr = mapBTables2model(BTable=BTable, model=model,
                                           allInter=F)
            newlinks = modelIntegr$reacID[modelIntegr$indexIntegr]
            """

    def run(self, model, data, k=2, compression=True,
            expansion=True, maxInputsPerGate=3, err1=0.1, err2=0,
            verbose=True):
        """

        :param model: filename to a SIF file
        :param data: filename to a MIDAS file
        :param int k: 2
        :param bool compression:
        :param bool expansion:
        :param int maxInputsPerGate:
        :param float err1:
        :param float err2:
        :param bool verbose:

        :return: nothing but populates the :attr:`newlinks` and
            :attr:`alllinks` attributes

        """
        from biokit.rtools import bool2R
        self.model = model
        self.data = data
        cmd = self.Rscript % {'pkn': model,
                              'data': data,
                              'k': k,
                              'compression': bool2R(compression),
                              'expansion': bool2R(expansion),
                              'err1': err1,
                              'err2': err2,
                              'maxInputsPerGate': maxInputsPerGate}

        if verbose:
            print(cmd)
        self.session.run(cmd)
        self.alllinks = self.session.get('alllinks')
        self.newlinks = self.session.get('newlinks')
        self.newlinks = [x for x in self.newlinks if x is not None]

    def __str__(self):
        txt = "Found {0} new links".format(len(self.newlinks))
        txt += str(self.newlinks)
        return txt

    def plot(self, cmap='jet', penwidth=5):
        """Plot graph with new links emphasized

        :param cmap: colormap
        :param penwidth: link width
        """
        from cno.io import CNOGraph
        from cno.io import Reactions

        c = CNOGraph(self.model, self.data)

        reactions = Reactions(self.newlinks)

        for e in c.edges():
            c.edge[e[0]][e[1]]['ecolor'] = 1

        for reaction in reactions._reactions:
            lhs, rhs = reaction.lhs, reaction.rhs
            # Could be a + or ^ combination
            for this_lhs in lhs.split("+"):
                if "^" not in this_lhs:
                    c.add_edge(this_lhs, rhs, link='+', ecolor=.5,
                               penwidth=3)  # all positive here by chance
                    c.edge[this_lhs][rhs]['penwidth'] = penwidth
                else:
                    # # TODO
                    raise NotImplementedError
        c.plot(edge_attribute='ecolor', edge_attribute_labels=False, cmap=cmap)
        return c
