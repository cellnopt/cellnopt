# -*- python -*-
#
#  This file is part of the cno package
#
#  Copyright (c) 2012-2014 - EMBL-EBI
#
#  File author(s): Thomas Cokelaer (cokelaer@ebi.ac.uk)
#
#  Distributed under the GLPv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: http://github.com/cellnopt/cellnopt
#
##############################################################################
from __future__ import print_function

import csv
import os
import re

from cno.io.reactions import Reactions, Reaction
from cno.misc import CNOError

import numpy as np

__all__ = ["SIF"]


class SIF(Reactions):
    """Manipulate network stored in SIF format.

    The SIF format is used in Cytoscape and CellNOpt (www.cellnopt.org).
    However, the format used in CellNOpt(R) restrict edges to be only
    1 or -1. Besides, special nodes called **AND** nodes can be added
    using the "and" string followed by a unique identifier(integer)
    e.g., and22; see below for details.

    .. seealso:: :ref:`sif` section in the online documentation.

    The SIF format is a tab-separated format. It encodes relations between nodes in a network.
    Each row contains a relation where the first column represents the input node, the second
    value is the type of relation. The following columns represents the output node(s). Here is a
    simple example::

        A 1 B
        B 1 C
        A -1 B

    but it can be factorised::

        A 1 B C
        B 1 C

    In SIF, only **OR** reactions can be encoded. The following::

        A 1 C
        B 1 C

    means A OR B gives C. **AND** reactions cannot be encoded therefore
    we have to code AND gates in a special way using a dedicated syntax.
    In order to encode the **AND** reaction the SIF reaction should be encoded
    as follows::

        A 1 and1
        B 1 and1
        and1 1 C

    An AND gate is made of the "and" string and a unique id concatenated as its end.

    A SIF file can be read as follows::

        s = SIF(filename)

    Each line is transformed into reactions (A=B, !A=B). You can then add or
    remove reactions. If you save the file in a new SIF file, be aware than
    lines such as::

        A 1 B C

    are expanded as::

        A 1 B
        A 1 C

    Aliases to the columns are stored in read-only attributes called :attr:`nodes1`,
    :attr:`edges`, :attr:`nodes2`. You can only add or remove reactions.
    Reactions are stored in :attr:`reactions`.



    """
    def __init__(self, filename=None, frmt="cno"):
        """.. rubric:: Constructor

        :param str filename: optional input SIF file.
        :param str frmt: "cno" or "generic" are accepted (default is cno). The cno format
            accepted only relation as "1" for activation, and "-1" for inhibitions. The "generic"
            format allows to have any relations. The "cno" format also interprets nodes that starts
            with "and" as logical AND gates. Such nodes are transformed into more cmpact notations.
            That is reactions.
        ::

            A 1 and1
            B -1 and1
            and1 1 C

        is transformed internally as 1 reaction: A^!B=C

        """
        super(SIF, self).__init__()

        self.frmt = frmt
        self.ignore_and = False
        if frmt == 'cno':
            self.convert_ands = True
        else:
            self.convert_ands = False

        # check that the extension is correct and the file exists
        if isinstance(filename, str):
            self.filename = filename
            if os.path.isfile(filename) is False:
                raise IOError("File %s not found." % filename)
            if filename.endswith('xml'):
                self.read_sbmlqual(filename)
            else:
                self.read_sif(filename)
        elif filename is None:
            self.clear()
        # could be another instance of a SIF file (or cnograph, or reactions)
        elif hasattr(filename, "reactions"):
            for reac in filename.reactions:
                self.add_reaction(reac)
        else:
            raise ValueError("argument must be a valid filename (string)")

    def clear(self):
        """remove all reactions and species"""
        self.remove_species(self.species)

    def is_and(self, value):
        if self.convert_ands is True and self.and_symbol in value:
            return True
        if re.search('^[a,A][n,N][D,d]\d+$', value) :
            return True
        return False

    def _get_and_nodes(self):
        _and_nodes = [x for x in set(self.nodes1 + self.nodes2)
                if self.is_and(x)]
        return _and_nodes
    and_nodes = property(_get_and_nodes, doc="Returns list of AND nodes")

    def remove_and_gates(self):
        """Remove all AND gates"""
        toremove = [r for x in self.and_nodes for r in self.reactions if x in r]
        for reac in toremove:
            self.remove_reaction(reac)

    def read_sif(self, filename):
        """Read a SIF file"""
        # Reads the data first
        self.clear()
        with open(filename, 'r') as fh:
            reader = csv.reader(fh)
            self._rawdata = [row for row in reader if len(row)]

        for i, row in enumerate(self._rawdata):
            if len(row[0].split()) < 3:
                raise Exception("Line %s contains a ill-formed reactions:%s" %(i+1, row))

        self._interpret_reactions()
        if self.convert_ands:
            self._convert_ands()

    def _convert_ands(self):
        #TODO check consistency between OR and AND gates

        # the and species found in the SIF (e.g A 1 and1; B 1 and1; and1 1 C)
        # are converted to A^B=C
        dont_remove = []
        for and_node in self.and_nodes:
            if self.and_symbol in and_node: # if ^ in name, nothing to do
                continue

            # otherwise, this is the original cno format that needs some mangling
            lhs_nodes = [(x,e) for x, e, y in zip(self.nodes1, self.edges, self.nodes2)
                    if y == and_node]
            rhs_node = [ y for x,e,y in zip(self.nodes1, self.edges,self.nodes2)
                    if x == and_node]

            try:
                assert len(rhs_node) == 1,  "%s %s %s" % (lhs_nodes, and_node, rhs_node)
                rhs_node = rhs_node[0]

                and_node = self.and_symbol.join([self.notedge(e)+x for x,e in lhs_nodes])
                reac =  and_node + "=" + rhs_node
                self.add_reaction(reac)
            except:
                dont_remove.append(and_node)

        # The and_node (e.g., and1) are finally removed
        for and_node in self.and_nodes:
            if and_node not in dont_remove and self.and_symbol not in and_node:
                self.remove_species(and_node)

    def _interpret_reactions(self):
        """interpret the data read from the SIF file"""
        for i, row in enumerate(self._rawdata):
            row = row[0].split() # row0 is the rawdata (string)
            node1 = row[0]
            edge = row[1]
            # SIF format allows several nodes on the RHS
            nodes = row[2:]
            for node2 in nodes:
                # some specific tests for CNO format
                if self.frmt == "cno":
                    if edge not in ["1","-1"]:
                        raise CNOError("Edges must be set to 1 or -1")
                    if re.search('^[a,A][n,N][D,d]\d+$',node1):
                        if edge == "-1":
                            raise ValueError("ill-formed SIF file line %s: an AND gate cannot have -1 edge" % i)

                # sometimes, we may want to ignore reactions with and reactions
                if self.ignore_and:
                    if self.is_and(node1) or self.is_and(node2):
                        continue

                # otherwise, store reactions
                if edge == "1":
                    reac = node1 + "=" + node2
                elif edge == "-1":
                    reac = "!" + node1 + "=" + node2
                else:
                    reac = node1 + "=" + node2
                self.add_reaction(reac)

    # makes the nodes1, nodes2, edges and data read-only properties
    def _get_nodes1(self):
        return [x.split("=")[0].replace("!","") for x in self.reactions]
    nodes1 = property(fget=_get_nodes1,
        doc="returns list of nodes in the left-hand sides of the reactions")

    def _get_nodes2(self):
        return [x.split("=")[1] for x in self.reactions]
    nodes2 = property(fget=_get_nodes2,
        doc="returns list of nodes in the right-hand sides of the reactions")

    def _get_edges(self):
        nodes1 = [x.split("=")[0] for x in self.reactions]
        edges =  ["-1" if x.startswith("!") else "1" for x in nodes1]
        return edges
    edges = property(fget=_get_edges,
        doc="returns list of edges found in the reactions")

    def add_reaction(self, reaction):
        """Adds a reaction into the network.

        Valid reactions are::

            A=B
            A+B=C
            A^B=C
            A&B=C

        Where the LHS can use as many species as desired. The following reaction
        is valid::

            A+B+C+D+E=F

        Note however that OR gates (+ sign) are splitted so A+B=C is added
        as 2 different reactions::

            A=C
            B=C

        """
        if "=" in reaction:
            lhs, rhs = reaction.split("=")
        else:
            raise ValueError("Reaction must contain a = sign")
        # OR gates can be splitted without issue. A+B=C can be added as 2
        # reactions A=C and B=C but AND gates are kept as it is with the symbol
        # ^
        if "+" in lhs:
            for specy in lhs.split("+"):
                reac = specy + "=" + rhs
                super(SIF, self).add_reaction(reac)
        else:
            reac = lhs + "=" + rhs
            super(SIF, self).add_reaction(reac)

    def save(self, filename):
        """Save the reactions (sorting with respect to order parameter)

        :param str filename: where to save the nodes1 edges node2

        """
        rhs = [x.rhs for x in self._reactions]
        # if we wer to use x.lhs, no need to select first item 
        # but ! are kept. So, we really need lhs_species.
        # It assumes there is only 1 item in the lhs, 
        # which should be true in the SIF format.
        lhs = [x.lhs_species[0] for x in self._reactions]
        
        f = open(filename, "w")
        sign2int = self.sign_operator_to_number

        counter = 1
        for reac in self._reactions:
            r = Reaction(reac)
            lhs_species = r.get_signed_lhs_species()
            if self.and_symbol in r.lhs:
                # create the and gate
                andname = "and%s" % counter
                for sign, species in lhs_species.items():
                    for this in species:
                        f.write("%s %s %s\n" % (this, sign2int(sign), andname))
                f.write("%s 1 %s\n" % (andname,r.rhs ))
                counter += 1
            else: # OR gate
                for sign, species in lhs_species.items():
                    for this in species:
                        f.write("%s %s %s\n" % (this, sign2int(sign), r.rhs))

        f.close()

    def sign_operator_to_number(self, operator):
        assert operator in ['+', '-']
        if operator == '+':
            return 1
        else:
            return -1

    def to_reactions(self):
        """Returns a Reactions instance generated from the SIF file.

        AND gates are interpreted. For instance the followinf SIF::

            A 1 and1
            B 1 and1
            and1 1 C

        give::

            A^B=C

        """
        from reactions import Reactions
        reactions = Reactions()
        for reac in self._reactions:
            reactions.add_reaction(reac)
        return reactions

    def notedge(self, x):
        """Returns ! character if x equals 1 and empty string otherwise"""
        if x=="-1":
            return "!"
        else:
            return ""

    def to_json(self):
        """Not a standard, do we want to keep this format ?

        """
        json = """{"links":[\n"""
        for n1, edge, n2 in zip(self.nodes1, self.edges, self.nodes2):
            json += """     {"source":%s, "target":%s, "link":%s},\n""" % (n1, n2, edge)
        json +="]}"
        return json

    def to_sbmlqual(self, filename=None):
        """Exports SIF to SBMLqual format.

        :param filename: save to the filename if provided
        :return: the SBML text

        This is a level3, version1 exporter.

        ::

            >>> s = SIF()
            >>> s.add_reaction("A=B")
            >>> res = s.to_sbmlqual("test.xml")

        .. warning:: logical AND are not encoded yet. works only if no AND gates

        .. warning:: experimental
        """

        #sif = self.to_cnograph()
        from cno.io.sbmlqual import SBMLQual
        qual = SBMLQual()
        sbml = qual.to_sbmlqual(self)

        if filename:
            fh = open(filename, "w")
            fh.write(sbml)
            fh.close()
        else:
            return sbml

    def read_sbmlqual(self, filename):
        """import SBMLQual XML file into a SIF instance

        :param str filename: the filename of the SBMLQual
        :param bool clear: remove all existing nodes and edges

        .. warning:: experimental
        """
        from cno.io.sbmlqual import SBMLQual
        qual = SBMLQual()
        sif = qual.read_sbmlqual(filename)
        self.clear()
        for reac in sif.reactions:
            self.add_reaction(reac)
        return sif

    def _rename_species(self, old, new):
        raise NotImplementedError
        for i, r in enumerate(self._reactions):
            self._reactions[i] = r.replace(old, new)

    def to_cnograph(self):
        # local import to prevent import cycling
        from cno.io.cnograph import CNOGraph
        c = CNOGraph()
        for reaction in self.reactions:
            c.add_reaction(reaction)
        return c

    def plot(self):
        """Plot the network


        .. note:: this method uses :class:`~cno.io.cnograph.CNOGraph` so
            AND gates appear as small circles.

        """
        c = self.to_cnograph()
        c.graph_options['graph']['dpi'] = 200
        c.plot()

    def __eq__(self, x):
        if isinstance(x, SIF) is False:
            return False
        if sorted(self.reactions) != sorted(x.reactions):
            return False
        return True

    def __str__(self):
        msg = "SIF object\n"
        msg += "- {0} reactions.\n".format(len(self.reactions))
        msg += "- {0} species.".format(len(self.species))
        msg += "\n"
        for reac in self.reactions:
            msg += reac + "\n"
        return msg

    def __len__(self):
        return len(self.nodes1)

    def __repr__(self):
        msg = "SIF object\n"
        msg += "- {0} reactions.\n".format(len(self.reactions))
        msg += "- {0} species.".format(len(self.species))
        return msg

