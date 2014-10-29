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

.. testsetup:: reactions

    from cno import Reaction
    from cno.io.reactions import Reactions
    r = Reactions()



"""
from __future__ import print_function
import re

import numpy
from easydev import Logging

from cno.misc import CNOError

__all__ = ["Reaction", "Reactions"]


class ReactionBase(object):
    valid_symbols = ["+","!", "&", "^"]


class Reaction(ReactionBase):
    """Logical Reaction

    A Reaction can encode logical ANDs and ORs as well as NOT::

        >>> from cno import Reaction
        >>> r = Reaction("A+B=C") # a OR reaction
        >>> r = Reaction("A^B=C") # an AND reaction
        >>> r = Reaction("A&B=C") # an AND reaction
        >>> r = Reaction("C=D")   # an activation
        >>> r = Reaction("!D=E")  # a NOT reaction

    The syntax is as follows:

    #. The **!** sign indicates a logical NOT.
    #. The **+** sign indicates a logical OR.
    #. The **=** sign indicates a relation (edge).
    #. The **^** or **&** signs indicate a logical AND. Note that **&** signs 
       will be replaced by **^**.


    Internally, reactions are checked for validity (e.g., !=C is invalid). 
    
    You can reset the name::

        >>> r.name = "A+B+C=D"

    or create an instance from another instance::

        >>> newr = Reaction(r)

    Sorting can be done inplace (default) or not. ::
    
        >>> r = Reaction("F+D^!B+!A=Z")
        >>> r.sort(inplace=False)
        '!A+!B^D+F=Z'

    Simple operator (e.g., equality) are available. The equality will sort the species
    internally so equality should be done on the instance (not the attribute :attr:`name`)::

        >>> r = Reaction("F+D^!B+!A=Z")
        >>> r.name == '!A+!B^D+F=Z'
        False
        >>> r == '!A+!B^D+F=Z'
        True

    If a reaction **A+A=B** is provided, it can be simplified by calling :meth:`simplify`.
    ANDs operator are not simplified. More sophisticated simplifications using Truth 
    Table could be used but will not be implemented in this class for now.
    """
    def __init__(self, reaction=None, strict_rules=True):
        """

        :param str reaction: a valid reaction (e.g., A=B, A+B=C, !B=D, C^D=F, ...),
            or an instance of :class:`Reaction`.
        :param bool strict_rules: if True, reactions cannot start with =, ^ or +
            signs (default to True).
        """
        self._strict_rules = strict_rules
        self._name = None
        if reaction is not None:
            # could be a Reaction instance
            if hasattr(reaction, "name"):
                self.name = reaction.name
            # or a string
            elif isinstance(reaction, str):
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

        .. doctest:: reactions

            >>> r = Reaction("!a+c^d^e^f+!h=b")
            >>> r.species
            ['a', 'c', 'd', 'e', 'f', 'h', 'b']
            

        """
        if reac == None:
            reac = self.name[:]
        species = re.split("[+|=|^|!]", reac)
        species = [x for x in species if x]
        return species
    species = property(_get_species)
    
    def _get_lhs(self):
        return self.name.split("=")[0]
    lhs = property(_get_lhs, 
            doc="Getter for the left hand side of the = character")

    def _get_lhs_species(self):
        lhs = self.name.split("=")[0]
        species = self._get_species(reac=lhs)
        return species
    lhs_species = property(_get_lhs_species,
            doc="Getter for the list of species on the left hand side of the = character")

    def _get_rhs(self):
        return self.name.split("=")[1]
    rhs = property(_get_rhs,
            doc="Getter for the right hand side of the = character")

    def _valid_reaction(self, reaction):
        reaction = reaction.strip()
        reaction = reaction.replace("&", "^")

        # = sign is compulsary
        N = reaction.count("=")
        if N != 1:
            raise CNOError("Invalid reaction name (only one = character expected. found %s)".format(N))
        #
        if self._strict_rules:
            if reaction[0] in ["=", "^" , "+"]:
                raise CNOError("Reaction (%s) cannot start with %s" %
                        (reaction, "=, ^, +"))

        # 
        lhs, rhs = reaction.split("=")
        for this in self.valid_symbols:
            if this in rhs:
                raise CNOError("Found an unexpected character (%s) in the LHS of reactions %s" % 
                        (reaction, self.valid_symbols))
        return reaction

    def sort(self, inplace=True):
        """Rearrange species in alphabetical order

        :param bool inplace: defaults to True

        ::

            >>> r = Reaction("F+D^!B+!A=Z")
            >>> r.sort()
            >>> r.name
            '!A+!B^D+F=Z'

        """
        # if only one lhs, nothing to do
        if len(self.lhs_species) == 1:
            return

        # we first need to split + and then ^
        splitted_ors = [x for x in self.lhs.split("+")] # left species keeping ! sign

        # loop over splitted list searching for ANDs
        species = []
        for this in splitted_ors:
            species_ands = this.split("^")
            # sort the species within the ANDs
            species_ands.sort(cmp=lambda x,y: cmp(x.replace("!", ""), y.replace("!", "")))
            species_ands = "^".join(species_ands)
            species.append(species_ands)

        # now sort the ORs
        species.sort(cmp=lambda x,y: cmp(x.replace("!", ""), y.replace("!", "")))
        # and finally rejoin them
        species = "+".join(species)

        new_reac = "=".join([species, self.rhs])
        if inplace:
            self.name = new_reac
        else:
            return new_reac

    def simplify(self, inplace=True):
        """Simplfies reaction if possible.

        ::

            >>> r = Reaction("A+A=B")
            >>> r.simplify()
            >>> r.name
            "A=B"

        Other cases (with ANDs) are not simplified.  Even though **A+A^B=C** truth table 
        could be simplified to **A=C** but we will not simplified it for now.

        """
        lhs = "+".join(set(self.lhs.split("+")))
        name = "=".join([lhs, self.rhs])
        if inplace:
            self.name = name
        else:
            return name

    def __eq__(self, other):
        # The reaction may not be sorted and user may not want to it to be sorted,
        # so we create a new instance and sort it
        r1 = Reaction(self)
        r1.sort()

        # we also sort the input reaction creating an instance as well so that the input reaction
        # (if it is an object) will not be sorted inplace either
        r2 = Reaction(other)
        r2.sort()
        if r1.name == r2.name:
            return True
        else:
            return False

    def __str__(self):
        txt = str(self.name)
        return txt


class Reactions(ReactionBase):
    """Data structure to handle list of :class:`Reaction` instances
    
    For the syntax of a reaction, see :class:`Reaction`. You can
    use the **=**, **!**, **+** and **^** characters.
    
    Reactions can be added using either string or instances of :class:`Reaction`::

        >>> from cno import Reaction, Reactions
        >>> r = Reactions()
        >>> r.add_reaction("A+B=C") # a OR reaction
        >>> r.add_reaction("A^B=C") # an AND reaction
        >>> r.add_reaction("A&B=C") # an AND reaction
        >>> r.add_reaction("C=D")   # an activation
        >>> r.add_reaction("!D=E")  # a NOT reaction
        >>> r.add_reaction(Reaction("F=G"))  # a NOT reaction

    Now, we can get the species::

        >>> r.species
        ['A', 'B', 'C', 'D', 'E']

    Remove one::

        >>> r.remove_species("A")
        >>> r.reactions
        ["B=C", "C=D", "!D=E"]

    .. note:: there is no simplifications made on reactions. For instance, if you add A=B and then 
        A+B=C, A=B is redundant but will be kept.

    .. seealso:: :class:`cno.io.reactions.Reaction` and :class:`cno.io.sif.SIF`
    """
    def __init__(self,  strict_rules=True, verbose=False):
        self._reactions = []
        self._reaction = Reaction(strict_rules=strict_rules)
        self.verbose = verbose
        self.strict_rules = strict_rules
       
    def to_list(self):
        """Return list of reaction names"""
        return [x.name for x in self.reactions]

    def _get_species(self):
        """Extract the specID out of reacID"""

        # extract species from all reacID and add to a set
        species = [this for reaction in self.reactions for this in reaction.species]
        species = set(species)

        # sort (transformed to a list)
        species = sorted(species)
        return species
    species = property(_get_species, doc="return list of unique species")

    def _get_reactions(self):
        return self._reactions
    reactions = property(fget=_get_reactions, doc="return list of reactions (objects)")

    def _get_reaction_names(self):
        return [reaction.name for reaction in self.reactions]
    reaction_names = property(fget=_get_reaction_names, doc="return list of reaction names")

    def __str__(self):
        _str = "Reactions() instance:\n" 
        _str += "- %s reactions\n" % len(self.reactions)
        _str += "- %s species\n" % len(self.species)
        return _str

    def remove_species(self, species_to_remove):
        """Removes species from the list of reactions

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
            lhs = reac.lhs_species  # lhs without ! sign
            rhs = reac.rhs

            # if RHS contains a species to remove, the entire reaction can be removed
            if rhs in species_to_remove:
                reacIDs_toremove.append(reac.name)
                continue

            # otherwise, we need to look at the LHS. If the LHS is of length 1,
            # we are in the first case (a=b) and it LHS contains specy to
            # remove, we do not want to keep it.
            if len(lhs) == 1:
                if lhs[0] in species_to_remove:
                    reacIDs_toremove.append(reac.name)
                    continue

            # Finally, if LHS contains 2 species or more, separated by + sign,
            # we do no want to remove the entire reaction but only the
            # relevant species. So to remove a in "a+b=c", we should return "b=c"
            # taking care of ! signs as well.
            for symbol in ["+", "^"]:
                if symbol not in reac.name:
                    continue
                else:
                    lhs_with_neg = [x for x in reac.name.split("=")[0].split(symbol)]
                    new_lhs = symbol.join([x for x in lhs_with_neg if x.replace("!", "") not in species_to_remove])
                    if len(new_lhs):
                        new_reac = new_lhs + "=" + rhs
                        reacIDs_toremove.append(reac.name)
                        reacIDs_toadd.append(new_reac)
        #
        for reac in reacIDs_toremove:
            self.remove_reaction(reac)

        for reac in reacIDs_toadd:
            self.add_reaction(reac)

    def add_reaction(self, reaction):
        """Adds a reaction in the list of reactions

        See documentation of the :class:`Reaction` for details. Here are
        some valid reactions::

            a=b
            a+c=d
            a^b=e  # same as above
            !a=e

        Example:
        
        .. doctest::

            >>> from cno import Reactions
            >>> c = Reactions()
            >>> c.add_reaction("a=b")
            >>> assert len(c.reactions) == 1

        """
        reac = Reaction(reaction, strict_rules=self.strict_rules)
        reac.sort()

        # 
        if reac.name not in self.to_list():
            self.reactions.append(reac)
        else:
            print("Reaction %s already in the list of reactions" % reaction)

    def remove_reaction(self, reaction_name):
        """Remove a reaction from the reacID list

            >>> c = Reactions()
            >>> c.add_reaction("a=b")
            >>> assert len(c.reactions) == 1
            >>> c.remove_reaction("a=b")
            >>> assert len(c.reactions) == 0
        """
        names = [x.name for x in self.reactions]
        if reaction_name in names:
            index2remove = names.index(reaction_name)
            del self.reactions[index2remove]
        else:
            if self.verbose:
                print("Reaction {0} not found. Nothing done".format(reaction_name))

    def search(self, species, strict=False):
        """Prints and returns reactions that contain the species name

        :param str species: name to look for
        :param bool strict: decompose reactions to search for the species
        :return: a Reactions instance with reactions containing the species to search for

        """
        r = Reactions()
        for x in self.reactions:
            list_species = x.lhs_species
            if strict == True:
                for this in list_species:
                    if species.lower() == this.lower():
                        if self.verbose:
                            print("Adding {0}".format(x.name))
                        r.add_reaction(x.name)
            else:
                for this in list_species:
                    if species.lower() in this.lower():
                        if self.verbose:
                            print("Adding {0}".format(x.name))
                        r.add_reaction(x.name)
                        continue
        return r

    def __len__(self):
        return len(self.reactions)

