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
"""This module provides tools to convert a SIF file into a format appropriate to
check sign consistency with ASP tools::

    A 1 B
    A -1 C

converted to ::

    A -> B +
    A -> C -

"""
from __future__ import print_function
#from __future__ import unicode_literals

from cno.io.sif import SIF

__all__ = ["SIF2ASP"]


class SIF2ASP(SIF):
    """Class to convert a SIF file into a ASP sign consistency format

    ::

        >>> from cno import SIF2ASP
        >>> from cno import cnodata
        >>> filename = cnodata("PKN-ToyMMB.sif")
        >>> s = SIF2ASP(filename)
        >>> s.to_net("PKN-ToyMMB.net")

    """
    def __init__(self, filename=None):
        """.. rubric:: Constructor

        :param str filename: the SIF filename

        """
        super(SIF2ASP, self).__init__(filename)

    def _get_signs(self):
        signs = ["+" if e=="1" else "-" for e in self.edges]
        return signs 
    signs = property(_get_signs, doc="get the signs of the reactions")

    def to_net(self, filename):
        """Write nodes and signs into a NET format

        If the SIF input format is ::

            A 1 B
            A -1 C

        the NET format should be::

            A -> B +
            A -> C -

        """
        if len(self.nodes1)>0:
            h = open(filename, "w")
            for n1,n2,s in zip(self.nodes1, self.nodes2, self.signs):
                h.write("%s -> %s %s\n" % (n1, n2, s))
            h.close()

    def __str__(self):
        txt = ""
        for n1,n2,s in zip(self.nodes1, self.nodes2, self.signs):
            txt += "%s -> %s %s" % (n1, n2, s)
        return txt



