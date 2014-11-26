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

from cno import CNOError

from cno.io.sbml import SBML
from cno.io.sif import SIF
from cno.io.cnograph import CNOGraph
from cno.io.reactions import Reaction, Reactions

import bs4

__all__ = ["SBMLQual"]


class SBMLQual(object):
    """Class to read and write SBML-qual file (logical models only)

    This is not an interface to SBML or SBML-qual. See libsbml library for
    that purpose. With this class you can read and write logical models
    stored in SBML-qual format.

    We do not guarantee that it covers all functionalities of SBML-qual but
    files saved with this class can be read back to be used within CellNOpt.

    You can convert CNOGraph of SIF instances to SBML-qual as follows::

    .. plot::
        :include-source:
        :width: 80%

        from cno import CNOGraph
        c1 = CNOGraph()
        c1.add_reaction("A+B=C")
        c1.expand_and_gates()
        c1.to_sbmlqual('test.xml')
        c1.plot()

        c2 = CNOGraph("test.xml')
        assert c1 == c2
        c2.plot()


    """
    def __init__(self):
        self.and_symbol = "^"

    def to_sbmlqual(self, graph):
        """Exports SIF to SBMLqual format.

        :return: the SBML text

        This is a level3, version 1 exporter.

        ::

            >>> s = SIF()
            >>> s.add_reaction("A=B")
            >>> res = s.to_SBMLQual("test.xml")

        """
        s = SBML(self, version="1.0", model_name='cellnopt_model')

        if isinstance(graph, SIF):
            data = graph.to_cnograph()
        elif isinstance(graph, CNOGraph):
            data = graph
        else:
            raise CNOError("Expected a CNOGraph of SIF intance as input")

        sbml = s.create_header()
        sbml += s.create_model_name()
        sbml += s.create_compartment(id="main", constant="true")

        # add the qualitativeSpecies list
        qualitativeSpecies = QualitativeSpecies(data)
        sbml += qualitativeSpecies.create()

        # Starting list of transitions
        list_of_transition = ListOfTransitions()
        sbml += list_of_transition.open()

        # Loop over all transitions
        tid = 0
        for node in sorted(data.nodes()):
            predecessors = data.predecessors(node)
            reactions = data.predecessors_as_reactions(node)

            if self.and_symbol in node:
                # if the output is a logical and, we skip the transition
                # it will be taken into account when it is an input
                continue
            if len(predecessors) == 0:
                continue # nothing to do, this is a source

            # else we have a new transition. We increment the identifier
            tid += 1
            identifier = "t{0}".format(tid)

            # and create a Transition
            transition = Transition(identifier)
            sbml += transition.open()

            # the list of inputs
            # - inputs could be an AND gate (e.g., A^!B=C), in which case, we want the
            #   species A and B to be extracted
            # - inputs could be made of positive and neg from same species (e.g., A=A
            #   and !A=A ), which means two entries in the list of inputs.
            species = {'-':[], '+':[]}
            for pred in predecessors:
                if self.and_symbol in pred:
                    d = Reaction(pred).get_signed_lhs_species()
                else:
                    if data[pred][node]['link'] == '+':
                        d = {'+': [pred]}
                    else:
                        d = {'-': [pred]}
                if '-' in d.keys():
                    species['-'].extend(d['-'])
                if '+' in d.keys():
                    species['+'].extend(d['+'])
            for k in species.keys():
                species[k] = list(set(species[k]))

            list_of_inputs = ListOfInputs(species, identifier)
            sbml += list_of_inputs.create()

            # The output (only one)
            list_of_outputs = ListOfOutputs(node)
            sbml += list_of_outputs.create()

            # Now the list of functions. This is the most complicated
            # but at the same time, we are lucky enough that in logical
            # models, the list of functions (at least in cellnopt) is
            # made of only one function, which is made of ORs. Inside
            # the ORs, you could have sveral ANDs
            list_of_function_terms = ListOfFunctionTerms()

            sbml += list_of_function_terms.open()
            sbml += list_of_function_terms.create_default_term()

            # there will be only one function term
            # if there is only one AND, starts with \and
            # else with \ors
            function_term = FunctionTerm()
            sbml += function_term.open()

            sbml +=  """<math xmlns="http://www.w3.org/1998/Math/MathML">"""
            if len(predecessors) == 1 and self.and_symbol in predecessors[0]:
                # a pure AND gate
                mathml = MathAND(predecessors[0], identifier)
                sbml += mathml.create()
            elif len(predecessors) == 1 and self.and_symbol not in predecessors[0]:
                # a direct link (no ORs, no ANDs)
                lhs = Reaction(reactions[0]).lhs
                sign = Reaction(reactions[0]).sign
                if sign == '1':
                    sign = '+'
                else:
                    sign = '-'
                    lhs = lhs.replace("!", "")
                mathml = MathApply(lhs, identifier, sign=sign)
                sbml += mathml.create()
            else: # an OR gate
                # inside the OR tag, you could have other gates
                # that is MathAND or MathApply
                # need to build a data structure that contains
                # the type of links. Needed for ORs only. ANDs
                # already contain the information in the name
                mathml = MathOR(reactions, identifier)
                sbml += mathml.create()
            sbml += "</math>"
            sbml += function_term.close()
            sbml += list_of_function_terms.close()
            sbml += transition.close()

        # The end
        sbml += list_of_transition.close()
        sbml += """</model>\n"""
        sbml += s.create_footer()

        sbml = self._prettify(sbml)

        return sbml

    def _prettify(self, sbml):
        """Return a pretty-printed XML string for the Element."""
        # beautifulsoup does a much better job than minidom but all tags are
        # transformed into lowercase.
        return bs4.BeautifulSoup(sbml).prettify()

        # not always the best layout
        #from xml.dom import minidom
        #reparsed = minidom.parseString(sbml)
        #return reparsed.toprettyxml(indent="    ", newl='')

    def read_sbmlqual(self, filename):
        """import SBMLQual XML file into a SIF instance

        :param str filename: the filename of the SBMLQual
        :param bool clear: remove all existing nodes and edges

        .. warning:: experimental
        """
        # We could just use XML.etree
        sif = SIF()
        res = bs4.BeautifulSoup(open(filename).read())

        # First, let us get the node names
        #model = res.findAll("model")[0]
        #allspecies = model.findAll("qual:listofqualitativespecies")[0]
        nodes  = [ x.get('qual:id') for x in res.findChildren("qual:qualitativespecies")]

        # Then, we go through all function terms
        for transition in res.findChildren("qual:transition"):
            inputs = [x['qual:qualitativespecies'] for x in transition.findChildren("qual:input")]
            signs = [x['qual:sign'] for x in transition.findChildren("qual:input")]
            output = [x['qual:qualitativespecies'] for x in transition.findChildren("qual:output")]
            assert len(output) == 1
            assert len(inputs) == len(signs)
            outputs = output * len(signs)

            # there may be different functions so we will need to loop over them
            functions = transition.findChildren("qual:functionterm")
            if len(functions)>1:
                CNOError("SBMLQual from cellnopt does not handle multiple functions")

            contents = functions[0].findChild('apply')
            if contents.find('and') and not contents.find('or'):
                lhs = self._get_lhs_from_apply(contents)
                reaction = lhs + "=" + outputs[0]
                sif.add_reaction(str(reaction))
            elif contents.find('or') and not contents.find('and'):
                lhs = self._get_lhs_from_apply(contents)
                reaction = lhs + "=" + outputs[0]
                sif.add_reaction(str(reaction))
            elif contents.find('or') is None and contents.find('and') is None:
                lhs = self._get_lhs_from_apply(contents)
                reaction = lhs + "=" + outputs[0]
                sif.add_reaction(str(reaction))
            else: #mulitple ORs
                for content in  contents.findChildren('apply', recursive=False):
                    lhs = self._get_lhs_from_apply(content)
                    reaction = lhs + "=" + outputs[0]
                    sif.add_reaction(str(reaction))

        # sanity check
        for node in nodes:
            if node not in sif.species:
                raise CNOError("A species without transition is not included in the network")

        return sif

    def _get_lhs_from_apply(self, xml):
        entries = xml.findChildren('apply', recursive=False)
        if len(entries) == 0:
            entries = [xml]
        lhs = []
        for entry in entries:
            if entry.find('geq') is not None:
                name = entry.find('ci').text.strip()
                lhs.append(name)
            else:
                name = entry.find('ci').text.strip()
                lhs.append("!" + name)
        if xml.find('and') is not None:
            lhs = "^".join(list(set(lhs)))
        else:
            lhs = "+".join(list(set(lhs)))
        return lhs


# NO NEED TO EXPORT ALL FOLLOWING CLASSES

# SBML-qual classes for logical modelling
class Qual(object):
    version = "http://www.sbml.org/sbml/level3/version1/qual/version1"

    def __init__(self, tag, xmlns=False):
        self.tag = tag
        self.xmlns = xmlns
        self.open_attribute = {}
        self.indent = ""

    def open(self):
        if self.xmlns is False:
            txt = """<qual:{0}""".format(self.tag)
            for k,v in self.open_attribute.items():
                txt+= """ qual:{0}="{1}" """.format(k,v) # note the space before 'qual'
            txt += ">\n"
        else:
            txt = """<qual:{0} xmlns:qual="{1}">""".format(self.tag, self.version)
        txt += "\n"
        return txt

    def close(self):
        return """</qual:{0}>\n""".format(self.tag)

    def indentation(self, sbml):
        sbml = "".join([self.indent + x for x in sbml.split("\n")])
        return sbml


class QualitativeSpecies(Qual):
    def __init__(self, species):
        super(QualitativeSpecies, self).__init__("listOfQualitativeSpecies", xmlns=True)
        self.species = species
        self.compartment = 'main'
        self.constant = 'false'

    def add_species(self, name):
        sbml = """<qual:qualitativeSpecies """
        sbml += """qual:constant="{0}" """.format(self.constant)
        sbml += """qual:compartment="{0}" """.format(self.compartment)
        sbml += """qual:id="{0}"/>\n""".format(name)
        return sbml

    def close(self):
        return ""

    def create(self):
        sbml = self.open()
        for name in self.species:
            if "^" not in name:
                sbml += self.add_species(name)
        sbml += self.close()
        sbml = self.indentation(sbml)
        return sbml


class ListOfTransitions(Qual):
    def __init__(self):
        super(ListOfTransitions, self).__init__("listOfTransitions", xmlns=True)


class Transition(Qual):
    """A transition contains at most one ListOfOnputs and one ListofOutputs and
    exactly one ListOfFunctionTerms

    A transition defines the level associated withthe QualitativeSpecies that occur
    when a Transition is enabled.

    In logical models a Transition is used to specify the logical rule associated with a
    QualitativeSpecies (that appears as an Output of this Transition). For example, the rule
    if A > 1: B = 2 would be encapsulated as a Transition with 2 QualitativeSpecies **A** as
    an input and **B**  as an Output; if A > 1 rule being encode by the math element of a
    3 FunctionTerm with the resultLevel attribute having a value 2.

    In Petri net models a Transition is interpreted, using the common Petri net
    semantics, as events that might occur within the system causing tokens to be moved.

    """
    def __init__(self, identifier):
        super(Transition, self).__init__("transition")
        self.identifier = identifier
        self.open_attribute = {'id':self.identifier}


class ListOfInputs(Qual):
    """The ListOfInputs contains at least one element of type Input.

    The input parameter **species** is a dictionay with keys + and - containing list
    of species in each category. A species could be in both categories.

    """
    def __init__(self, species, identifier):
        super(ListOfInputs, self).__init__("listOfInputs")
        self.species = species
        self.identifier = identifier
        assert '+' in self.species.keys()
        assert '-' in self.species.keys()
        self.threshold = 1
        self.transitionEffect = 'none'

    def create(self):
        txt = self.open()

        # positive and then negative:
        prefix = """<qual:input qual:thresholdLevel="{0}" """.format(self.threshold)
        prefix += """ qual:transitionEffect="{0}" """.format(self.transitionEffect)
        for name in self.species['+']:
            txt += prefix
            txt += """ qual:sign="positive" """
            txt += """ qual:qualitativeSpecies="{0}" """.format(name)
            txt += """ qual:id="theta_{0}_{1}"/>""".format(self.identifier, name)
        for name in self.species['-']:
            txt += prefix
            txt += """ qual:sign="negative" """
            txt += """ qual:qualitativeSpecies="{0}" """.format(name)
            txt += """ qual:id="theta_{0}_{1}"/>""".format(self.identifier, name)

        txt += self.close()
        return txt


class ListOfOutputs(Qual):
    """In logical model, there is only one output

    * thresholdLevel is set to 1
    * transitionEffect is set to assignmentLevel


    """
    def __init__(self, node):
        super(ListOfOutputs, self).__init__('listOfOutputs')
        self.name = node

    def create(self):
        txt = self.open()
        txt += """<qual:output qual:thresholdLevel="1" """
        txt += """ qual:transitionEffect="assignmentLevel" """
        txt += """ qual:qualitativeSpecies="{0}"/>\n""".format(self.name)
        txt += self.close()
        return txt


class ListOfFunctionTerms(Qual):
    """

    contains 1 default terms and any number of function terms

    """
    def __init__(self):
        super(ListOfFunctionTerms, self).__init__('listOfFunctionTerms')

    def create_default_term(self):
        default = DefaultTerm()
        return default.create()

    #def add_list_function_term(self):
    #    raise NotImplementedError


class FunctionTerm(Qual):
    """associated with a result and to a boolean function inside a math element
    that can be used to set the conditions inder which this term is selected


    """
    def __init__(self):
        super(FunctionTerm, self).__init__('functionTerm')
        self.open_attribute = {'resultLevel': '1'}


class MathApply(object):
    def __init__(self, name, identifier, sign="+"):
        self.name = name
        self.identifier = identifier
        assert sign in ['+', '-']
        self.sign = sign

    def create(self):
        txt = "<apply>\n"
        if self.sign == '+':
            txt += "<geq/>\n"
        else:
            txt += "<lt/>\n"
        txt += "<ci> {0} </ci>\n".format(self.name)
        txt += "<ci> theta_{0}_{1} </ci>\n".format(self.identifier, self.name)
        txt += "</apply>\n"
        return txt


class MathOR(object):
    def __init__(self, reactions, identifier):
        self.reactions = Reactions(reactions)
        self.identifier = identifier

    def create(self):
        txt = '<apply>\n'
        txt += '<or/>\n'
        for reaction in self.reactions._reactions:
            if "^" in reaction.name:
                ml = MathAND(reaction.name, self.identifier)
                txt += ml.create()
            else:
                if reaction.sign == '1':
                    sign = '+'
                else:
                    sign = '-'
                name = reaction.lhs
                name = name.replace("!", "")
                ml = MathApply(name, self.identifier, sign)
                txt += ml.create()
        txt += '</apply>\n'
        return txt


class MathAND(object):
    """Get MathML representation of an AND gate.

    """
    def __init__(self, reaction, identifier):
        """
        identifier is the transition identifier
        """
        self.reaction = Reaction(reaction)
        self.identifier = identifier

    def create(self):
        txt = '<apply>\n'
        txt += '<and/>\n'
        species = self.reaction.get_signed_lhs_species()
        for name in species['+']:
            mathapply = MathApply(name, self.identifier)
            txt += mathapply.create()
        for name in species['-']:
            mathapply = MathApply(name, self.identifier, '-')
            txt += mathapply.create()
        txt += '</apply>\n'
        return txt


#class Math(Qual):
#    def __init__(self):
#        super(Math,self).__init__('math')##
#
#    def open(self):
#        return """<math xmlns="http://www.w3.org/1998/Math/MathML">"""


class DefaultTerm(Qual):
    """resultLevel is set to 0"""
    def __init__(self):
        super(DefaultTerm, self).__init__('defaultTerm')
        self.open_attribute = {'resultLevel': 0}

    def create(self):
        txt = self.open()
        txt += self.close()
        return txt

