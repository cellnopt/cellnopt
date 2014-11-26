from cno.io.reactions import Reaction
import pandas as pd
import itertools

__all__ = ['TruthTable']


class TruthTable(object):
    def __init__(self, reaction):

        if "=" not in reaction:
            reaction += "="
        self.reaction = Reaction(reaction) 
        self.reaction.sort() # sorted name is in self.reaction.name 
        # TODO inside reaction, sort inplace once created

        self.species = sorted(set(self.reaction.lhs_species))
        #assert len(set(self.species)) == len(self.species)

        self.ors = self.reaction.lhs.split("+")
        
        table = list(itertools.product([1, 0], 
            repeat=len(self.species)))

        # Create the dataframe for the species first so we
        # have all the permutation.
        self.df = pd.DataFrame(table, columns=self.species)

        # adds each individual part that is complex (^ sign)
        for this in self.ors:
            if "^" in this:
                self._eval_ands(this)
            else:
                self._eval_species(this)

        # compute the final reaction
        values = self.df[self.reaction.lhs.split("+")].max(axis=1)
        self.df[self.reaction.lhs] = values
        self.df.columns = sorted(list(self.df.columns), key=lambda x: len(x))

    def _get_tt(self):
        return self.df[self.species+[self.reaction.lhs]]
    truthtable = property(_get_tt)

    def __eq__(self, other):
        # species are sorted to build the table
        # so this shjoudl suffice
        # TODO: check carefully
        return all(self.truthtable == other.truthtable)

    def _eval_species(self, species):
        """Evaluate A or !A if A present in self.df"""
        assert species in self.species or species[1:]

        if species not in self.df.columns:
            if species.startswith("!"):
                single_name = species[1:]
                self.df[species] = 1. - self.df[single_name]

    def _eval_ands(self, reaction):
        assert "^" in reaction
        assert "+" not in reaction
        # make sure all species are known
        for species in reaction.split("^"):
            print("    " + species)
            self._eval_species(species)

        # now compute the ANDs.
        columns = reaction.split("^")
        self.df[reaction] = list(self.df[columns].min(axis=1))





