# -*- python -*-
#
#  This file is part of cellnopt software
#
#  Copyright (c) 2014-2014 - EBI-EMBL
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
"""ASP related"""
import pandas as pd
from cno.io.cnograph import CNOGraph
from cno.io.sif import SIF

__all__ = ["NET", "net2reaction"]


class NET(object):
    """Class to manipulate reactions in NET format.

    The NET format ::

        species1 -> species2 sign

    where sign can be either the + or - character.


    Examples are::

            A -> B +
            A -> B -



    """
    def __init__(self, filename=None):
        """.. rubric:: constructor

        :param str filename: optional filename containing NET reactions
            if provided, NET reactions are converted into reactions (see
            :class:`cellnopt.core.reactions.Reactions`

        """
        self._net = []
        self.read_net(filename=filename)

    def read_net(self, filename=None):
        """read NET file


        """
        if filename is None:
            return

        try:
            f = open(filename, 'r')
            data = f.read()
            data = data.splitlines()
            for i,row in enumerate(data):
                if len(row)>0:
                    self.add_net(row)
                else:
                    print("warning. found an empty line. skipped")
        except IOError as err:
            raise IOError(err)
        else:
            # if the file was opened, let us close it.
            f.close()

    def _get_reactions(self):
        sif = self.to_sif()
        return sif.reactions
    reactions = property(_get_reactions)

    def to_sif(self, filename=None):
        """Write SIF reactions into a file"""
        sif = SIF()
        for net in self.net:
            sif.add_reaction(self.net2reaction(net))
        if filename:
            sif.save(filename)
        else:
            return sif

    def _get_net(self):
        return self._net
    net = property(_get_net)

    def add_net(self, net):
        self._check(net)
        self._net.append(net)

    def _check(self, net):
        try:
            self._split_net(net)
        except Exception as err:
            raise Exception(err)

    def __str__(self):
        txt = ""
        for row in self.net:
            txt += row + "\n"
        return txt

    def _split_net(self, data):
        try:
            lhs, rhs = data.split("->")
        except:
            raise ValueError("net2reaction: it seems your net string is ",
                             "not correct. could not find -> characters")
        lhs = lhs.strip()
        try:
            rhs, sign = rhs.split()
            rhs = rhs.strip()
            sign = sign.strip()
        except:
            raise ValueError("RHS of your net string could not be split. missing  space ? ")
        if sign not in ['+', '-']:
            raise ValueError("found invalid sign. Must be either + or - character")
        return (lhs, rhs, sign)

    def net2reaction(self, data):
        """convert a NET string to a reaction

        a NET string can be one of ::

            A -> B +
            C -> D -

        where + indicates activation and - indicates inhibition

        .. doctest::

            >>> assert net2reaction("A -> B +") == "A=B"
            >>> assert net2reaction("A -> B -") == "!A=B"

        """
        lhs, rhs ,sign = self._split_net(data)
        if sign == "+":
            reaction = "=".join([lhs, rhs])
        elif sign == "-":
            reaction = "!" + "=".join([lhs, rhs])
        return reaction


"""
class CASPOModels(object):

    def __init__(self, filename):
        print("Deprecated. USe Models class instead.")
        self.filename = filename
        self.df = pd.read_csv(self.filename)
        self.df.columns = [x.replace("+", "^") for x in self.df.columns]

        self.cnograph = CNOGraph()
        for this in self.df.columns:
            self.cnograph.add_reaction(this)

    def get_average_model(self):
        return self.df.mean(axis=0)


    def plotdot(self, model_number=None, *args, **kargs):
        if model_number==None:
            model = self.get_average_model()
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
            self.cnograph.edge[edge[0]][edge[1]]["label"] = value
            self.cnograph.edge[edge[0]][edge[1]]["average"] = value

        self.cnograph.plotdot(edge_attribute="average", **kargs)

    def to_sif(self, filename):

        self.cnograph.to_sif(filename)
"""

