from cno.misc.truthtable import TruthTable

def test_truthtable():


    # could a reaction with strict rules or not (rhs present or not)
    t1 = TruthTable("A+B")
    t2 = TruthTable("A+B=C")
    assert t1 == t2

    t1 = TruthTable("A^!B")
    t2 = TruthTable("!B^A=C")
    assert t1 == t2
