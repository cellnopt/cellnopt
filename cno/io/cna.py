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
""":Topic: **Module dedicated to the CNA reactions data structure**


:Status: for production but not all features implemented.

"""
from __future__ import print_function


from cno.io.reactions import Reactions
from cno.io.sif import SIF
from cno.misc import CNOError

__all__ = ["CNA"]


class CNA(Reactions):
    """Reads a reaction file (CNA format)

    This class has the :class:`Interaction` class as a Base class.
    It is used to read **reactions** files from the CNA format, which 
    is a CSV-like format where each line looks like::

        mek=erk   1 mek = 1 erk   |   #  0 1 0   436  825  1    1  0.01 

    The pipe decompose the strings into a LHS and RHS. 

    The LHS is made of a unique identifier without blanks (mek=erk). The remaining part is
    the reaction equation. The equal sign "=" denotes the reaction arrow. Identifiers, 
    coefficients and equal sign must be separated by at least one blank. 
    The ! sign to indicate not. The + sign indicates an OR relation. 

    .. warning:: The + sign indicates an OR as it should be. However, keep in
        mind that in CellNOptR code, the + sign
        indicates an AND gate. In this package we always use **+** for an OR and 
        **^** or **&** for an AND gate. 

    .. warning:: in the CNA case, some reactions have no LHS or RHS. Such
        reactions are valid in CNA but may cause issue if converted to SIF
        

    .. note:: there don't seem to be any AND in CNA reactions.

    The RHS is made of

    * a default value: # or a value. 
    * a set of 3 flags representing the time scale 
       * flag 1: whether this interaction is to be excluded in logical computations
       * flag 2: whether the logical interaction is treated with incomplete truth table
       * flag 3: whether the interaction is monotone
    * reacBoxes (columns 5,6,7,8)
    * monotony  (col 9)

    In this class, only the LHS are used for now, however, the RHS values are
    stored in different attributes.

    ::

        >>> from cno.io import Reactions
        >>> from cno import getdata
        >>> a = Reactions(getdata('test_reactions'))
        >>> reacs = a.reactions

    .. seealso:: CNA class inherits from :class:`cno.io.reaction.Reaction`
    """
    def __init__(self, filename=None, type=2, verbose=False):
        """.. rubric:: Constructor


        :param str filename: an optional filename containing reactions in CNA
            format. If not provided, the CNA object is empty but you can
            add reactions using :meth:`~cno.io.cna.CNA.add_reaction`. 
            However, attributes such as :attr:`~cno.io.cna.CNA.reacBoxes` 
            will not be populated.
        :param integer type: only type 2 for now. 
        :param bool verbose: False by default

        .. todo:: type1 will be implemented on request.

        """
        super(CNA, self).__init__()
        self.strict_rules = False

        #self.metabolites = metabolites # must be a class LoadMetabolites
        self.filename = filename
        self.verbose = verbose
        self.type = type
        if type != 2:
            raise NotImplementedError("only type 2 implemented")


        # Attributes populated while reading the data.
        #: populated when reading CNA reactions file
        self.reacBoxes = []
        #: populated when reading CNA reactions file
        self.incTruthTable = []
        #: populated when reading CNA reactions file
        self.timeScale = []
        #: populated when reading CNA reactions file
        self.excludeInLogical = []
        #: populated when reading CNA reactions file
        self.reacText = []
        #: populated when reading CNA reactions file
        self.monotony = []  #flag 3
        self.reacDefault = []
 
        if filename:
            self._read_reactions()
            self._get_species()

    def _read_reactions(self):
        """Read a reactions file and populate readID"""
        f = open(self.filename, "r")
        data = []                   # the data structure to populate
        for line in f.readlines():  # for each line
            # convert tab.to white space, remove trailing and \n character
            line = line.replace('\t',' ').replace('\n','').strip()

            # do not consider commented or empty lines
            if line.startswith("%") or line.startswith('#'): 
                pass
            if len(line) == 0:
                print("Found an empty line. Skipped")
            else:
                data.append(line)
        f.close()

        # scan all the data
        for i, x in enumerate(data):
            try:
                beforePipe, afterPipe = x.split('|') # there should be only one pipe per
                # line, so if it fails, this is a format error
            except ValueError as err:
                raise ValueError("Error msg to do")

            reacID = beforePipe.split()[0].strip()
            if reacID.count('=') != 1:
                raise ValueError("Error line %s: wrong format expected one " %(i+1)
                    + "only one = sign, found %s" % reacID.count('='))
            else:
                self.add_reaction(reacID)
            reacText = beforePipe.replace(reacID, "").strip()
            self.reacText.append(reacText)
            parameters = afterPipe.split()
            if len(parameters) != 9:
                raise ValueError("Error line %s: did no find expected numbers of parameters" % i+1)

            if self.type == 1:
                # not finished
                reacDefault, reacMin, reacMax, objFunc, d, d, d, d, reacVariance = parameters
                mue = []
                stoichMat = []
            elif self.type == 2:
                # First, the reac default value.
                if parameters[0].isalnum():
                    self.reacDefault.append(float(parameters[0]))
                elif parameters[0].strip()=='#':
                    self.reacDefault.append(float('NaN'))
                else:
                    raise ValueError("""Error line %s: unexpected value in the
first column after pipe character (%s)""" % (str(i+1), parameters[0]))

                self.incTruthTable.append(float(parameters[1]))
                self.timeScale.append(float(parameters[2]))
                self.excludeInLogical.append(float(parameters[3]))
                self.monotony.append(float(parameters[8]))
                self.reacBoxes.append([i+1, float(parameters[4]), 
                    float(parameters[5]), 0, float(parameters[6]),
                    float(parameters[7])])
            # clean up the reacDefault: could be # or number
        if self.verbose == True:
            print(self)

    def to_sif(self, filename=None):
        """Export the reactions to SIF format 

        ::

            from cno.io import CNA
            r = CNA()
            r.add_reaction("a=b")
            r.add_reaction("a+c=e")
            r.to_sif("test.sif")

        Again, be aware that "+" sign in Reaction means "OR". 
        Looking into the save file, we have the a+c=e reactions (a=e OR c=e)
        expanded into 2 reactions (a 1 e) and (c 1 e) as expected::

            a   1   b
            a   1   e
            c   1   e

        """
        s = SIF()
        for reac in self.reactions:
            try:
                s.add_reaction(reac)
            except CNOError:
                print("Skipped {} reaction".format(reac))
        s.save(filename)
