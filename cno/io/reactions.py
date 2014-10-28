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
"""This module contains a base class to manipulate reactions

"""
from __future__ import print_function

from easydev import Logging

import numpy

from cno.misc import CNOError

__all__ = ["Reaction"]



class ReactionBase(object):
    valid_symbols = ["+","!", "&", "^"]


class Reaction(ReactionBase):
    """A Reaction class 

    A Reaction can encode logical AND and OR as well as NOT::

        >>> from cellnopt.core import Reaction
        >>> r = Reaction("A+B=C") # a OR reaction
        >>> r = Reaction("A^B=C") # an AND reaction
        >>> r = Reaction("A&B=C") # an AND reaction
        >>> r = Reaction("C=D")   # an activation
        >>> r = Reaction("!D=E")  # a NOT reaction

            r.name
            r.rename_species(old, new)
            r._valid_reaction("a=b")

    .. todo:: what is the meaning of A+A=B. ERROR, WARNING, SILENT by removing duplicates?
        A+!A=B ? A^!A=B
    """
    def __init__(self, reaction=None, strict_rules=True):
        """

        :param str reaction:
        :param bool strict_rules: if True, reactions cannot start with =, ^ or ^
            signs.
        """
        self._strict_rules = strict_rules
        self._name = None
        if reaction is not None:
            self.name = reaction

    def _set_name(self, reaction):
        if reaction is not None:
            reaction = self._valid_reaction(reaction)
        self._name = reaction[:]
    def _get_name(self):
        return self._name
    name = property(_get_name, _set_name, doc="Getter/Setter for the reaction name")

    def _get_species(self, reac=None):
        """

            >>> r = Reaction("!a+c^d^e^f+!h=b")
            >>> r._get_species()
            ['a' ,'c', 'd', 'r' ,'f' ,'h' ,'b']
        """
        if reac==None:
            reac = self.name[:]
        import re
        species = re.split("[+|=|^|!]", reac)
        species = [x for x in species if x]
        return species
    species = property(_get_species)
    
    def _get_lhs(self):
        return self.name.split("=")[0]
    lhs = property(_get_lhs)

    def _get_lhs_species(self):
        lhs = self.name.split("=")[0]
        species = self._get_species(reac=lhs)
        return species
    lhs_species = property(_get_lhs_species)

    def _get_rhs(self):
        return self.name.split("=")[1]
    rhs = property(_get_rhs)

    def _valid_reaction(self, reaction):
        reaction = reaction.strip()
        reaction = reaction.replace("&", "^")

        #
        if self._strict_rules:
            if reaction[0] in ["=", "^" , "+"]:
                raise CNOError("Reaction (%s) cannot start with %s" %
                        (reaction, "=, ^, +"))
        # = sign is compulsary
        if "=" not in reaction:
            raise CNOError("Reaction (%s) must contain a = sign" % reaction)

        # 
        lhs, rhs = reaction.split("=")
        for this in self.valid_symbols:
            if this in rhs:
                raise CNOError("Found an unexpected character (%s) in the LHS of reactions %s" % (reaction, self.valid_symbols))
        if "&" in lhs or "^" in lhs:
            if lhs == "":
                raise CNOError("Found an AND gate without RHS (%s)" % (reaction))
        return reaction

    def __str__(self):
        if self.name: 
            return self.name
        else: 
            return ""



class Reactions(ReactionBase):
    """A class to manipulate list of reactions (e.g., A=B)
    
    You can create list of reactions using the **=**, **!**, **+** 
    and **^** characters with the following meaning::

        >>> from cellnopt.core import *
        >>> c = Reactions()
        >>> c.add_reaction("A+B=C") # a OR reaction
        >>> c.add_reaction("A^B=C") # an AND reaction
        >>> c.add_reaction("A&B=C") # an AND reaction
        >>> c.add_reaction("C=D")   # an activation
        >>> c.add_reaction("!D=E")  # a NOT reaction

     #. The **!** sign indicates a logical NOT.
     #. The **+** sign indicates a logical OR.
     #. The **=** sign indicates a relation (edge).
     #. The **^** or **&** signs indicate a logical AND  **&** signs will be 
        replaced by **^**.

    .. warning:: the **+** sign being a logical OR, A+B=C is the same as adding 
       two reactions: A=C, B=C. 

    Now, we can get the species::

        >>> c.species
        ['A', 'B', 'C', 'D', 'E']

    Remove one::

        >>> c.remove_species("A")
        >>> c.reactions
        ["B=C", "C=D", "!D=E"]

    .. seealso:: :class:`cellnopt.core.reactions.Reactions` and :class:`cellnopt.core.sif.SIF`
    """

    def __init__(self, format="cno", strict_rules=True):
        self._reactions = []
        self._format = format
        self._logging = Logging("INFO")
        self._reaction = Reaction(strict_rules=strict_rules)
       
    def to_list(self):
        return [x.name for x in self.reactions]

    def _reac2spec(self, reac):
        """simple function to extract species from a reaction """
        self._reaction._valid_reaction(reac)
        reac = reac.replace("!","") # we don't care about the NOT here
        reac = reac.replace("^","+")  # just to find the species
        lhs, rhs = reac.split("=")

        specID = []
        specID.extend([x for x in lhs.split('+') if x])
        specID.extend([x for x in rhs.split('+') if x])
        specID = set(specID)
        return specID

    def _get_species(self):
        """Extract the specID out of reacID"""
        specID = set() # use set to prevent duplicates

        # extract species from all reacID and add to the set
        for reac in self.reactions:
            specID = specID.union(self._reac2spec(reac))

        # sort (transformed to a list)
        specID = sorted(specID)
        return specID
    species = property(_get_species, doc="return species")

    def _get_reactions(self):
        return self._reactions
    reactions = property(fget=_get_reactions)

    def _get_lhs_species(self, reac, remove_not=True):
        """Return independent species

        reac = "A+!B+C=D"
        _get_lhs_species returns ['A', 'B', 'C']

        used by _build_interMat
        """
        lhs = reac.split('=')[0]  # keep only lhs
        lhs = lhs.split('+')      # split reactions
        lhs = [x for x in lhs if len(x)>0]
        if remove_not == True:
            lhs = [x.replace('!','') for x in lhs]
        return lhs

    def __str__(self):
        _str = "Reactions() instance:\n" 
        _str += "- %s reactions\n" % len(self.reactions)
        _str += "- %s species\n" % len(self.species)
        return _str

    def remove_species(self, species_to_remove):
        """Removes species from the reactions list

        :param str,list species_to_remove:


        .. note:: If a reaction is "a+b=c" and you remove specy "a", 
            then the reaction is not enterely removed but replace by "b=c"

        """
        # make sure we have a **list** of species to remove
        if isinstance(species_to_remove, list):
            pass
        elif isinstance(species_to_remove, str):
            species_to_remove = [species_to_remove]
        else:
            raise TypeError("species_to_remove must be a list or string")

        reacIDs_toremove = []
        reacIDs_toadd = []
        for reac in self.reactions:
            lhs = self._get_lhs_species(reac)  # lhs without ! sign
            rhs = reac.split("=")[1]

            # two cases: either a=b or a+d+e=b
            # if we remove a, the first reaction should be removed but the
            # second should just be transformed to d+e=b

            # RHS contains a specy to remove, we do want the reaction
            if rhs in species_to_remove:
                reacIDs_toremove.append(reac)
                continue

            # otherwise, we need to look at the LHS. If the LHS is of length 1,
            # we are in the first case (a=b) and it LHS contains specy to
            # remove, we do not want to keep it.
            if len(lhs) == 1:
                if lhs[0] in species_to_remove:
                    reacIDs_toremove.append(reac)
                    continue

            # Finally, if LHS contains 2 species or more, separated by + sign,
            # we do no want to remove the entire reaction but only the
            # relevant specy. So to remove a in "a+b=c", we should return "b=c"
            # taking care of ! signs.
            for symbol in ["+", "^"]:
                if symbol not in reac:
                    continue
                else:
                    lhs_with_neg = [x for x in reac.split("=")[0].split(symbol)]
                    new_lhs = symbol.join([x for x in lhs_with_neg if x.replace("!", "") not in species_to_remove])
                    if len(new_lhs):
                        new_reac = new_lhs + "=" + rhs
                        reacIDs_toremove.append(reac)
                        reacIDs_toadd.append(new_reac)

        level = self._logging.level
        self._logging.level = "ERROR"
        for reac in reacIDs_toremove:
            self.remove_reaction(reac)
        for reac in reacIDs_toadd:
            self.add_reaction(reac)
        self._logging.level = level


    def add_reaction(self, reaction):
        """Adds a reaction in the list of reactions

        In logical formalism, the inverted hat stand for OR but there is no such
        key on standard keyboard so we use the + sign instead. The AND is
        defined with either the ^ or & sign. Finally the NOT is 
        defined by the ! sign. Valid reactions are therefore::

            a=b
            a+c=d
            a&b=e
            a^b=e  # same as above
            !a=e

        Example::

            >>> c = Reactions()
            >>> c.add_reaction("a=b")
            >>> assert len(c.reactions) == 1

        .. warning & symbol are replaced by ^ internally.

        """
        reac = Reaction(reaction)
        #reaction = self._sort_reaction(reaction)

        # TODO 
        if reac.name not in self.to_list():
            self.reactions.append(reac)
        else:
            self._logging.info("%s already in the list of reactions" % reaction)

    def _sort_reaction(self, reaction):
        """Rearrange species of the LHS in alphabetical order

        ::

            >>> s.sort_reaction("b+a=c")
            "a+b=c"

        """
        lhs, rhs = reaction.split("=")
        lhs = [x for x in lhs.split("+")] # left species with ! sign   ["b","!a","c"]
        if len(lhs) == 1:
            return reaction
        # if more than one specy, we must rearrange them taking care of ! signs
        species = self._get_lhs_species(reaction)  # without ! signs ["b", "a", "c"]

        sorted_indices = numpy.argsort(species)
        new_lhs = []
        for i in sorted_indices:
            new_lhs.append(lhs[i])
        new_reac = "=".join(["+".join(new_lhs), rhs])
        return new_reac

    def remove_reaction(self, reaction):
        """Remove a reaction from the reacID list

            >>> c = Reactions()
            >>> c.add_reaction("a=b")
            >>> assert len(c.reactions) == 1
            >>> c.remove_reaction("a=b")
            >>> assert len(c.reactions) == 0
        """
        if reaction in self.reactions:
            self._reactions.remove(reaction)

    def search(self, species, strict=False, verbose=True):
        """Prints and returns reactions that contain the species name

        :param str species: name to look for
        :param bool strict: decompose reactions to search for the species
        :return: a Reactions instance with reactions containing the species to search for

        """
        r = Reactions()
        for x in self.reactions:
            species = self._reac2spec(x)
            print(species)
            if strict == True:
                for this in species:
                    if species.lower() == this.lower():
                        if verbose:
                            print(x)
                        r.add_reaction(x)
            else:
                for this in species:
                    if species.lower() in this.lower():
                        if verbose:
                            print(x)
                        r.add_reaction(x)
                        continue
        return r



