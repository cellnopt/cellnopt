# -*- python -*-
#
#  This file is part of the cinapps.tcell package
#
#  Copyright (c) 2012-2014 - EMBL-EBI
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

import csv
import os
import re

# could be replace since we just use 2 functions ?
from numpy import sort, array

from cno.io.reactions import Reactions

__all__ = ["SIF"]


class SIF(Reactions):
    """Manipulate network stored in SIF format.

    The SIF format is used in Cytoscape and CellNOpt (www.cellnop.org).
    However, the format used in CellNOpt(R) restrict edges to be only
    1 or -1. Besides, special nodes called **AND** nodes can be added
    using the "and" string followed by a unique identifier(integer)
    e.g., and22; seebelow for details.

    .. seealso:: :ref:`sif` section in the online documentation.

    The SIF format is a tab-separated format. It encodes relations betwee nodes in a network.
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
    Reactions are stored in :attr:`reacID`.


    .. todo:: explain more precisely or simplify the 2 parameter ignore_and
        and convert_ands, which are different semantic ! one of r the ^ character,
        one for the and string.

    """
    def __init__(self, filename=None, format="cno", ignore_and=False, convert_ands=True):
        """.. rubric:: Constructor

        :param str filename: optional input SIF file.
        :param str format: "cno" or "generic" are accepted (default is cno). The cno format
            accepted only relation as "1" for activation, and "-1" for inhibitions. The "generic"
            format allows to have any relations. The "cno" format also interprets nodes that starts
            with "and" as logical AND gates.
        :param bool ignore_and: if you want to ignore the and nodes (see above), set to True.
        :param bool convert_ands: if AND nodes are found (from cellnopt syntax, eg a^b), converts them
            into a single reaction (default is True).
        """
        super(SIF, self).__init__()

        self.format = format
        self.ignore_and = ignore_and
        self.convert_ands = convert_ands
        self.clear()

        # check that the extension is correct and the file exists
        if isinstance(filename, str):
            self.filename = filename
            #extension = os.path.splitext(filename)[1]
            #if extension not in [".sif", ".SIF"]:
            #    raise ValueError("SIF file must have extension sif or SIF. %s found." % extension)
            if os.path.isfile(filename) == False:
                raise IOError("File %s not found." % filename)
            self.loadSIF(filename)
        elif filename==None:
            pass
        elif hasattr(filename, "reactions"):
            for reac in filename.reactions:
                self.add_reaction(reac.name)

        else:
            raise ValueError("argument must be a valid filename (string)")

    #def clear(self):
    #    # simply reset the relevant attributes
    #    self._edges = []
    #    self._nodes1 = []
    #    self._nodes2 = []
    #    self._data = []
    #    self._rawdata = []
    #    for reac in self.reacID:
    #        self.remove_reaction(reac)
    def clear(self):
        """remove all reactions and species"""
        self.remove_species(self.species)

    def _get_and_nodes(self):
        _andNodes = [x for x in set(self.nodes1 + self.nodes2) if
                re.search('^[a,A][n,N][D,d]\d$',x)]
        if self.convert_ands:
            _andNodes += [x for x in set(self.nodes1 + self.nodes2) if "^" in x]
        return _andNodes
    andNodes = property(_get_and_nodes, doc="Returns list of AND nodes")

    def remove_and_gates(self):
        toremove = [r for x in self.andNodes for r in self.reaction_names if x in r]
        for reac in toremove:
            self.remove_reaction(reac)

    def loadSIF(self, filename):
        # Reads the data first
        self.clear()
        try:
            f = open(filename, 'r')
            reader = csv.reader(f)
            self._rawdata = [row for row in reader if len(row)]

            for i, row in enumerate(self._rawdata):
                if len(row[0].split()) < 3:
                    raise Exception("Line %s contains a ill-formed reactions:%s" %(i+1, row))
        except IOError, e:
            raise IOError(e)
        else:
            # if the file was opened, let us close it.
            f.close()

        self._interpret_reactions()
        self._convert_ands()

    def _convert_ands(self):
        #TODO check consistency between OR and AND gates
        if self.convert_ands == False:
            return

        # the and species found in the SIF (e.g A 1 and1; B 1 and1; and1 1 C)
        # are converted to A^B=C
        dont_remove = []
        for andNode in self.andNodes:
            if "^" in andNode:
                continue
            lhs_nodes = [ (x,e) for x,e,y in zip(self.nodes1, self.edges,self.nodes2) if y==andNode]
            rhs_node = [ y for x,e,y in zip(self.nodes1, self.edges,self.nodes2) if x==andNode]
            try:
                assert len(rhs_node) == 1,  "%s %s %s" % (lhs_nodes, andNode, rhs_node)
                rhs_node = rhs_node[0]

                andNode = "^".join([self.notedge(e)+x for x,e in lhs_nodes])
                reac =  andNode + "=" + rhs_node
                self.add_reaction(reac)
            except:
                dont_remove.append(andNode)

        # The andNode (e.g., and1) are finally removed
        for andNode in self.andNodes:
            if andNode not in dont_remove and "^" not in andNode:
                self.remove_species(andNode)

    def _interpret_reactions(self):
        """interpret the data read from the SIF file"""
        for i, row in enumerate(self._rawdata):
            row = row[0].split() # row0 is the rawdata (string)
            node1 = row[0]
            edge = row[1]
            nodes = row[2:] # may be several nodes after edge (at least 1)
            for node2 in nodes:
                # some specific tests for CNO format
                if self.format == "cno":
                    if edge not in ["1","-1"]:
                        raise ValueError("Edges must be set to 1 or -1")
                    if re.search('^[a,A][n,N][D,d]\d+$',node1):
                        if edge == "-1":
                            raise ValueError("ill-formed SIF file line %s: an AND gate cannot have -1 edge" % i)

                # sometimes, we may want to ignore reactions with and reactions
                if self.ignore_and:
                    if "and" in node1 or "and" in node2:
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
        return [x.split("=")[0].replace("!","") for x in self.reaction_names]
    nodes1 = property(fget=_get_nodes1,
        doc="returns list of nodes in the left-hand sides of the reactions")

    def _get_nodes2(self):
        return [x.split("=")[1] for x in self.reaction_names]
    nodes2 = property(fget=_get_nodes2,
        doc="returns list of nodes in the right-hand sides of the reactions")

    def _get_edges(self):
        nodes1 = [x.split("=")[0] for x in self.reaction_names]
        edges =  ["-1" if x.startswith("!") else "1" for x in nodes1]
        return edges
    edges = property(fget=_get_edges,
        doc="returns list of edges found in the reactions")

    def _get_data(self):
        data = array([(x,y,z) for x,y,z in
            zip(self.nodes1, self.edges, self.nodes2)],
            dtype=[('nodes1',object),('edges','int'),('nodes2',object)])
        # instead of object, we can use a string type S but we must provide the
        # length. So, we need to estimate the lnogest specy name first.
        return data
    data = property(fget=_get_data, doc="Returns list of relations")

    def __eq__(self, x):
        if isinstance(x, SIF) == False:
            return False
        if sorted(self.nodes1) != sorted(x.nodes1):
            return False
        if sorted(self.nodes2) != sorted(x.nodes2):
            return False
        if sorted(self.edges) != sorted(x.edges):
            return False
        return True

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

    def save(self, filename, order="nodes1"):
        """Save the reactions (sorting with respect to order parameter)

        :param str filename: where to save the nodes1 edges node2
        :param str order: can be nodes1, edges or nodes2

        """
        assert order in ["nodes1", "nodes2", "edges"]
        f = open(filename, "w")
        for row in sort(self.data, order=order):
            f.write("%s %s %s\n" % (row[0], str(row[1]), row[2]))
        f.close()

    def __str__(self):
        """prints in alphabetical order"""
        msg = self._underline("Original file") + "\n\n"
        M = max([len(x) for x in self.nodes1+self.nodes2])
        M2 = max([len(x) for x in self.edges])
        for row in sort(self.data, order="nodes1"):
            msg += "%*s %*s %*s\n" % (M+2, row[0], M2+2 , str(row[1]), M+2,row[2])
        msg = msg.rstrip("\n")

        if self.format == "cno":
            msg += "\n\n"
            msg += self._underline("namesSpecies")+"\n"
            msg += ", ".join([x for x in self.species]) + "\n\n"
            msg += self._underline("reacID")+"\n"
            msg += ", ".join([x for x in self.reaction_names]) + "\n\n"
        return msg

    def _underline(self, msg):
        import easydev
        return easydev.underline(msg)


    def sif2reaction(self):
        """Returns a Reactions instance generated from the SIF file.

        AND gates are interpreted. For instance the followinf SIF::

            A 1 and1
            B 1 and1
            and1 1 C

        give::

            A^B=C

        """
        reacID = []
        for n1, e, n2 in zip(self.nodes1, self.edges, self.nodes2):
            if n1 not in self.andNodes and n2 not in self.andNodes:
                if e == "-1":
                    reacID.append("!" + "=".join([n1,n2]))
                else:
                    reacID.append("=".join([n1,n2]))

        for andNode in self.andNodes: # what to do with and nodes ?
            lhs_nodes = [ (x,e) for x,e,y in zip(self.nodes1, self.edges,self.nodes2) if y==andNode]
            rhs_node = [ y for x,e,y in zip(self.nodes1, self.edges,self.nodes2) if x==andNode]
            assert len(rhs_node) == 1, rhs_node
            rhs_node = rhs_node[0]

            reac = "^".join([self.notedge(e)+x for x,e in lhs_nodes])
            reac +=  "="+rhs_node
            reacID.append(reac)

        from reactions import Reactions
        reactions = Reactions()
        for reac in reacID:
            reactions.add_reaction(reac)
        return reactions

    def notedge(self, x):
        """Returns ! character if x equals 1 and empty string otherwise"""
        if x=="-1":
            return "!"
        else:
            return ""

    def sif2json(self):
        raise NotImplementedError

    def to_SBMLQual(self, filename=None, modelName="CellNOpt_model"):
        """Exports SIF to SBMLqual format.

        :param filename: save to the filename if provided
        :param str modelName: name of the model in SBML document
        :return: the SBML text

        This is a level3, version1 exporter.

        ::

            >>> s = SIF()
            >>> s.add_reaction("A=B")
            >>> res = s.to_SBMLQual("test.xml")

        .. warning:: logical AND are not encoded yet. works only if no AND gates

        .. warning:: experimental
        """
        from cellnopt.core import SBML
        s = SBML(self, version="1.0", model_name=modelName)

        if len(self.andNodes)>=1:
            print("and gates are handled safely. Use with great care !!!")
            print("!!!! READ THE COMMENT ABOVE")

        sbml = s.create_header()
        sbml += s.create_model_name()
        sbml += s.create_compartment(id="main", constant="true")


        # add the species first
        sbml += """\n    <qual:listOfQualitativeSpecies xmlns:qual="http://www.sbml.org/sbml/level3/version1/qual/version1">"""
        for node in self.species:
             sbml += """\n      <qual:qualitativeSpecies qual:constant="false" qual:compartment="main" qual:id="%s"/>""" % node
        sbml += """\n    </qual:listOfQualitativeSpecies>\n"""

        c = self.to_cnograph()

        sbml += """\n    <qual:listOfTransitions xmlns:qual="http://www.sbml.org/sbml/level3/version1/qual/version1">"""
        tid = 0
        for node in sorted(c.nodes()):
            if "^" in node: continue  # otherwise pickup in AND and OR cases...
            predecessors = c.predecessors(node)
            if len(predecessors):
                tid +=1
                transitionName = "t%s" % tid

                # the inputs
                sbml += """\n      <qual:transition qual:id="%s">""" % transitionName
                sbml += """\n        <qual:listOfInputs>"""
                for pred in predecessors:
                    if "^" not in pred:
                        sign = c[pred][node]['link']
                        if sign == "+":  sign = "positive"
                        elif sign == "-":  sign = "negative"
                        else: raise ValueError("sign (%s) if wrong expected + or - sign)" % sign)
                        sbml += """
          <qual:input qual:thresholdLevel="1" qual:transitionEffect="none" qual:sign="%s" qual:qualitativeSpecies="%s" qual:id="theta_%s_%s"/>""" % (sign, pred, transitionName, pred)
                    else:
                        input_ands = pred.split('=')[0].split("^")
                        for this in input_ands:
                            if this.startswith("!"):
                                pred = pred[1:]
                                sign = 'negative'
                            else:
                                sign= 'positive'
                            sbml += """
          <qual:input qual:thresholdLevel="1" qual:transitionEffect="none" qual:sign="%s" qual:qualitativeSpecies="%s" qual:id="theta_%s_%s"/>""" % (sign,this, transitionName, this)

                sbml += """\n        </qual:listOfInputs>\n"""

                # the output (only one)
                sbml += """\n        <qual:listOfOutputs>"""
                sbml += """\n          <qual:output qual:transitionEffect="assignmentLevel" qual:qualitativeSpecies="%s"/>""" % node
                sbml += """\n        </qual:listOfOutputs>\n"""

                sbml +="""\n        <qual:listOfFunctionTerms>
          <qual:defaultTerm qual:resultLevel="0">
          </qual:defaultTerm>
          <qual:functionTerm qual:resultLevel="1">
            <math xmlns="http://www.w3.org/1998/Math/MathML">"""

                # the equations depending on the number of inputs and type of
                # edge (-1 or 1)
                if len(predecessors)==1 and "^" not in predecessors[0]:
                    link = c[predecessors[0]][node]['link']
                    if link == "+": relation = "<geq/>"
                    elif link == "-": relation = "<lt/>"
                    else: raise ValueError("wrong sign must be + or -")
                    sbml +="""
              <apply>
                %s
                <ci> %s </ci>
                <ci> theta_%s_%s </ci>
              </apply>""" % (relation, predecessors[0], transitionName, predecessors[0])
                elif len(predecessors)==1 and "^" in predecessors[0]:
                    input_ands = predecessors[0].split('=')[0].split("^")
                    sbml += """
            <apply>
              <and/>"""
                    for pred in input_ands:
                        if pred.startswith("!"):
                            pred = pred[1:]
                            relation = "<lt/>"
                        else:
                            relation = "<geq/>"
                        sbml +="""
                <apply>
                  %s
                  <ci> %s </ci>
                  <ci> theta_%s_%s </ci>
                </apply>""" % (relation, pred, transitionName, pred)

                    sbml += """\n              </apply>"""


                elif len(predecessors)>1:
                    for pred in predecessors:
                        if "^" in pred:
                            raise ValueError("combination of inputs must be made of simple ORs without ANDS. case not yet implemented")
                    sbml += """
              <apply>
                <or/>"""
                    for pred in predecessors:
                        link = c[pred][node]['link']
                        if link == "+": relation = "<geq/>"
                        elif link == "-": relation = "<lt/>"
                        else: raise ValueError("wrong sign must be + or -")
                        sbml +="""
                <apply>
                  %s
                  <ci> %s </ci>
                  <ci> theta_%s_%s </ci>
                </apply>""" % (relation, pred, transitionName, pred)

                    sbml += """\n              </apply>"""
                else:
                    print("Unknown case")
                sbml += """\n            </math>"""


            if len(predecessors)>=1:
                sbml += """\n          </qual:functionTerm>"""
                sbml += """\n        </qual:listOfFunctionTerms>"""
                sbml += """\n      </qual:transition>\n"""


        # The end
        sbml += """    </qual:listOfTransitions>\n"""
        sbml += """  </model>\n"""
        sbml += s.create_footer()

        if filename:
            fh = open(filename, "w")
            fh.write(sbml)
            fh.close()

        return sbml

    def importSBMLQual(self, filename, clear=True):
        """import SBMLQual XML file into a SIF instance

        :param str filename: the filename of the SBMLQual
        :param bool clear: remove all existing nodes and edges

        .. warning:: experimental
        """
        if clear:
            self.clear()
        import bs4
        # !!!!!!!!!!!!!!!!!!!!!!! somehow xml is transformed in lower case?????
        res = bs4.BeautifulSoup(open(filename).read())
        model = res.findAll("model")[0]
        allspecies = model.findAll("qual:listofqualitativespecies")[0]

        nodes  = [ x.get('qual:id') for x in res.findChildren("qual:qualitativespecies")]

        andCounter = 1
        for transition in res.findChildren("qual:transition"):
            inputs = [x['qual:qualitativespecies'] for x in transition.findChildren("qual:input")]
            signs = [x['qual:sign'] for x in transition.findChildren("qual:input")]
            output = [x['qual:qualitativespecies'] for x in transition.findChildren("qual:output")]
            assert len(output) == 1
            assert len(inputs) == len(signs)
            outputs = output * len(signs)

            # there may be different functions so we will need to loop over them
            functions = transition.findChildren("qual:functionterm")

            for function in functions:
                # from function, figure out if we have a OR or AND gate
                ands = function.findChildren("and")
                #TODO naive approach. could use somthing more robust ???
                # should include the OR and other function terms that are not AND
                if len(ands)==0:
                    print(inputs, signs, outputs)
                    for x,link,y in zip(inputs, signs, outputs):
                        if link == "positive":
                            self.add_reaction("%s=%s" % (x,y))
                        else:
                            self.add_reaction("!%s=%s" % (x,y))

                # should include the AND. there is no else because I',m not sure if
                # there could both OR and AND equations withi a functionTerms.
                if len(ands):
                    print("Found an AND gate")
                    andName = "and%s" % andCounter
                    LHS = []
                    for x,link in zip(inputs, signs):
                        if link == "positive":
                            LHS.append("%s" % x)
                        else:
                            LHS.append("!%s" % x)

                    reaction = "^".join(LHS) + "=%s" % outputs[0]
                    print(reaction)
                    self.add_reaction(reaction)
                andCounter += 1

        for node in nodes:
            if node not in self.species:
                raise NotImplementedError("A species without transition is not included in the network")

    def _rename_species(self, old, new):
        print("Draft method. Use with care.")
        # could be mced to  Interactions class
        for i, r in enumerate(self._reactions):
            self._reactions[i] = r.replace(old, new)

    def to_cnograph(self):
        from cellnopt.core import CNOGraph
        c = CNOGraph()
        for reaction in self.reaction_names:
            c.add_reaction(reaction)
        return c

    def plot(self):
        """Plot the network


        .. note:: this method uses :class:`cellnopt.core.cnograph` so AND gates appear
            as small circles.

        """
        c = self.to_cnograph()
        c.graph_options['graph']['dpi'] = 200
        c.plot()

    def __len__(self):
        return len(self.nodes1)
