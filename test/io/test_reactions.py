from cno.io import Reaction, Reactions


def test_reaction():
    r = Reaction("!a+c^d^e^f+!h=b")
    r.lhs_species
    r.lhs
    r.rhs
    print(r)

    r = Reaction("A=B")
    try:
        r = Reaction("=B")
        assert False
    except:
        assert True
    
    # this is possible with the strict_rule set on
    r = Reaction("=B", strict_rules=False)

    try:
        r = Reaction("A")
        assert False
    except:
        assert True

    try:
        r = Reaction(1)
        assert False
    except:
        assert True

    try:
        r = Reaction("A++B=C")
        assert False
    except:
        assert True

    try:
        r = Reaction("A!B=C") # missing + or ^ before ! sign
        assert False
    except:
        assert True


    # simplification
    r = Reaction("A+A=B")
    r.simplify()
    assert r.name == "A=B"

    # simplification
    r = Reaction("A+A=B")
    assert r.simplify(inplace=False) == "A=B"

    # sorting
    r = Reaction("F+D^!B+!A=Z")
    r.sort()
    assert r.name == '!A+!B^D+F=Z'
    r = Reaction("F+D^!B+!A=Z")
    assert r.sort(inplace=False) == '!A+!B^D+F=Z'

    assert (Reaction("A=B") == "C=B") is False

    assert Reaction("A=B") != Reaction("C=B")


    r = Reaction("A+B^!C=D")
    assert r.name == "A+B^!C=D"
    r.rename_species({'A':'F', 'C': 'E', 'D':'DD'})
    assert r.name == 'F+B^!E=DD'

    # test get_signed 
    r = Reaction("A+B^!C=D")
    r.get_signed_lhs_species()
    r = Reaction("A+!B=D")
    r.get_signed_lhs_species()
    r = Reaction("A+B=D")
    r.get_signed_lhs_species()


    r1 = Reaction("A+B=C")
    r2 = Reaction("B+A=C")
    assert r1 == r2



class test_Reactions():

    def _get_reactions(self):
        return [Reaction('A=B'), Reaction('C=B')]

    def test_rename_species(self):
        rs = Reactions(self._get_reactions())
        rs.rename_species({'A':'a', 'B':'b'})
        assert sorted(rs.reactions) == ['C=b', 'a=b']


    def test_reactions(self):
        c = Reactions(verbose=True)
        assert c.reactions == []
        assert c.species == []
        c.add_reaction("A=B")
        c.add_reaction("A=B")
        assert len(c) == 1
        c.reactions
        c.remove_reaction("A=B")
        assert len(c) == 0
        c.remove_reaction("A=B")

    def test_wrong_reactions(self):
        c = Reactions()
        try:
            c.add_reaction("A")
            assert False
        except:
            assert True
        for symbol in c.valid_symbols:
            try:
                c.add_reaction("A=%sB" % symbol)
                assert False
            except:
                assert True        

    def test_search(self):
        c = Reactions(verbose=True)
        c.add_reaction("MEK=ERK")
        c.search("MEK")
        c.search("ME", strict=True)
        c.search("MEK", strict=True)


    def test_remove_species(self):
        c = Reactions()
        c.add_reaction("A+B=C")
        c.remove_species("A")
        c.reactions == ['B=C']


        c.add_reaction("A=C")
        c.remove_reaction("A=C")
        c.reactions == ['B=C']

        # no effect
        c.remove_species("dummy")

        # should be empty after this
        c.remove_species(["A", "B", "C"])
        assert len(c.reactions) == 0

        c.add_reaction("A=B")
        c.remove_species("A")

        # wrong type
        try:
            c.remove_species(1)
            assert False
        except:
            assert True



