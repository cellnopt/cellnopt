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
from __future__ import print_function

from  cno.io.reactions import Reactions


__all__ = ["SOP2SIF"]


class SOP2SIF(Reactions):
    r"""Converts a file from SOP to SIF format

    SOP stands for sum of products, it is a list of relations of the form::

        !A+B=C

    For now, this function has been tested and used on the 
    copy/paste of a PDF document
    into a file. Be careful because the interpretation of the characters may differ 
    from one distribution to the other. The original data contains

        #. a special character for NOT, which is interpreted as \x2\xac (a L turned by 90 degrees clockwise)
        #. an inversed ^ character for OR, which is interpreted as " _ "
        #. a ^ character for AND, which is correctly interpreted.
        #. a -> character for "gives", which is transformed into ! character.

    On other systems, it may be interpreted differently, so we provide a mapping
    attribute :attr:`mapping` to perform the translation, which can be changed to your needs.

    The data looks like::

        1 !A + B = C 1 [references]
        2 !A + B = E 2 [references]
        3 !A + B = D 1 [references]
        ...
        N !A + B = D 2 [references]

    The :class:`SOP2SIF` class gets rid of the last column, the [references] and the column
    before it (made of 1 and 2). Then, we convert the reaction strings into the
    same format as in CellNOpt that is:

        #. A = C means A GIVES C
        #. A + B = C means A gives C OR B gives C 
        #. !A   means NOT A

    ::

        >>> s2s = SOP2SIF("data.sop")
        >>> s = s2s.sop2sif()
        >>> s2s.writeSIF("data.sif")


    """
    def __init__(self, filename):
        super(SOP2SIF, self).__init__()
        self.filename = filename
        #self.data = None
        self._read_data()
        self.or_string = "__or__"
        # values of the mapping directory must be changed to match the contents
        # of the filename provided.

        #: The dictionary to map SOP special characters e.g if you code NOT with ! character, just fill this dictionary accordingly
        self.mapping = {
            'not':'\xc2\xac',
            'gives':'!',
            'and':'^',
            'or':' _ '}
        self._translate()

    def _read_data(self):
        """Reads the data and performs some cleanup"""

        # reads the data
        fh = open(self.filename, "r")
        self.data = fh.read()
        fh.close()

        # split lines
        self.data = self.data.split("\n")[0:-1]

        # strip spaces
        self.data = [x.strip() for x in self.data]

        # remove line with comments and empty lines
        self.data = [x for x in self.data if not x.startswith('#') and len(x)>0]

    def _translate(self):
        """Interprets the data and translate into proper format.

        This function transform all lines into proper reactions expanding the
        AND and OR gates.

        A^B=C (AND gate) is transformed into 3 reactions::

            A=and1
            B=and1
            and1=C

        The OR gates (e.g., A+B=C) are split as well::

            A=C
            B=C

        """
        if self.mapping['not'] != '\xc2\xac' \
                or self.mapping['gives'] != '!' \
                or self.mapping['or'] != ' _ ':
            raise NotImplementedError

        data = self.data[:]

        # replace the direction character into equal. This should be done
        # before the NOT because in some cases, the direction sign is coded as
        # !, which is the NOT sign. 
        data = [x.replace(self.mapping['gives'], '=') for x in data]
        data = [x.replace('= ', '=') for x in data]
        data = [x.replace(' =', '=') for x in data]

        # replace a special character by ! and get rid of spaces
        #if self.mapping['='] == '!':
        data = [x.replace(self.mapping['not'], '!') for x in data]
        data = [x.replace('! ', '!') for x in data]
        data = [x.replace(' ! ', '!') for x in data]

        # replace the "and" character by + and remove spaces if any
        data = [x.replace(self.mapping['and'], '^') for x in data]
        data = [x.replace(' ^', '^') for x in data]
        data = [x.replace('^ ', '^') for x in data]
     
        # get rid of parantheses
        data = [x.replace('(', '') for x in data]
        data = [x.replace(')', '') for x in data]

        # remove spaces, and replace the or by a special tag to be interpreted
        # later by the writeSIF method.
        if self.mapping['or'] == " _ ":
            #data = [x.replace(" _ ", self.or_string) for x in data]
            data = [x.replace(" _ ", "+") for x in data]
        else:
            raise NotImplementedError


        # Remove the first column (line number)
        data = [x.split(" ", 1)[1] for x in data]

        # If the reactions has no spaces, we can get rid of the
        # the tau and reference columns by simply performing a split operation.
        data = [x.split(" ")[0] for x in data]

        # some negative reactions are written as A=!B, which is
        # identical to !A=B, which is our convention.
        for i, reac in enumerate(data):
            lhs, rhs = reac.split('=')
            if rhs.startswith('!'):
                lhs = '!'+lhs
                rhs = rhs.replace('!','')
                newreac = lhs+'='+rhs
                print("Warning: found a ! in RHS. inversion performed: %s is now %s" %(rhs, newreac))
                data[i] = newreac

        # cleanup the or gates by splitting them: A__or__B=C becomes A=C and B=C
        # and start to fill reactions
        ORreacs = [x for x in data if "+" in x]
        data = [x for x in data if "+" not in x]

        # add back the OR reactions 
        for reacs in ORreacs:
            lhs, rhs = reacs.split("=")
            for l in lhs.split("+"):
                newreac = l + '=' + rhs
                data.append(newreac)

        # cleanup the AND gates by creating AND nodes
        ANDreacs = [x for x in data if "^" in x]
        data = [x for x in data if "^" not in x]
        for i,reac in enumerate(ANDreacs):
            andNode = "and%s" % str(i+1)
            lhs, rhs = reac.split("=")
            species = lhs.split("^")
            for specy in species:
                data.append("%s=%s" % (specy, andNode))
            data.append("%s=%s" % (andNode, rhs))



        self._reactions = []
        for r in data:
            lhs, rhs = r.split('=')
            if "^" in rhs or "+" in rhs:
                print("Warning: reaction %s skipped (several + or ^ signs in RHS)" % r),
                continue
            self.add_reaction(r)

        print("Parsing done. Found %s reactions." % len(self.species))


    def sop2sif(self, include_and_gates=True):
        """Converts the SOP data into a SIF class


        :param bool include_and_gates: if set to False, 
            all reactions with AND gates are removed.

        :returns: an instance of :class:`cno.io.sif.SIF`

        """
        from cno.io.sif import SIF
        s = SIF()
        for reac in self.reactions:
            s.add_reaction(reac)
        if include_and_gates == False:
            s.remove_and_gates()
        return s

    def export2sif(self, filename, include_and_gates=True):
        """Save the reactions in a file using SIF format

        The data read from the SOP file is transformed into a SIF class before
        hand.

        :param bool include_and_gates: if set to False, all reactions with AND
            gates removed
        """
        s = self.sop2sif(include_and_gates)
        s.save(filename)


    def __len__(self):
        return len(self.reactions)
