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
from cno.io.sbml import SBML


__all__ = ["SBMLQual"]


class SBMLQual(object):
    def __init__(self):
        pass

    def to_sbmlqual(self, sif):
        """Exports SIF to SBMLqual format.

        :return: the SBML text

        This is a level3, version1 exporter.

        ::

            >>> s = SIF()
            >>> s.add_reaction("A=B")
            >>> res = s.to_SBMLQual("test.xml")

        .. warning:: logical AND are not encoded yet. works only if no AND gates

        .. warning:: experimental
        """
        # FIXME use bs4 to make the code nicer.

        s = SBML(self, version="1.0", model_name='cellnopt_model')

        if len(sif.andNodes)>=1:
            print("and gates may not be handled safely. Use with great care !!!")
            print("!!!! READ THE COMMENT ABOVE")

        sbml = s.create_header()
        sbml += s.create_model_name()
        sbml += s.create_compartment(id="main", constant="true")

        # add the species first
        sbml += """\n    <qual:listOfQualitativeSpecies xmlns:qual="http://www.sbml.org/sbml/level3/version1/qual/version1">"""
        for node in sif.species:
             sbml += """\n      <qual:qualitativeSpecies qual:constant="false" qual:compartment="main" qual:id="%s"/>""" % node
        sbml += """\n    </qual:listOfQualitativeSpecies>\n"""

        c = sif.to_cnograph()

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

        return sbml

    def read_sbmlqual(self, filename):
        """import SBMLQual XML file into a SIF instance

        :param str filename: the filename of the SBMLQual
        :param bool clear: remove all existing nodes and edges

        .. warning:: experimental
        """
        from cno.io.sif import SIF
        sif = SIF()
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
                    for x,link,y in zip(inputs, signs, outputs):
                        if link == "positive":
                            sif.add_reaction("%s=%s" % (x,y))
                        else:
                            sif.add_reaction("!%s=%s" % (x,y))

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
                    sif.add_reaction(reaction)
                andCounter += 1

        for node in nodes:
            if node not in sif.species:
                raise NotImplementedError("A species without transition is not included in the network")

        return sif

