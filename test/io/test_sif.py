from cno.io.sif import SIF
from cno.testing import getdata
import easydev
import tempfile

simple = getdata("test_simple.sif")
double = getdata("test_double.sif")
net = getdata("PKN-ToyMMB.net")
  


def test_sifreader():
    """Test the 2 different ways of coding SIF files

    double contains a relation:
     
        A 1 B C

    that is coded in simple as:

        A 1 B
        A 1 C

    """
    s1 = SIF(simple)
    s2 = SIF(double)

    assert s1.nodes1 == ['A', 'B', 'A', 'E']
    assert s1.nodes2 == ['B', 'C', 'C', 'F']
    assert s1.edges == ['1', '1', '1', '-1']

    assert s1 == s2

    assert s1.species == s2.species
    print(s1)
    len(s1)

    # we can also create an instance from another instance:
    s1 = SIF(s2)
    assert  s1 == s2

    # let us call the clear method:
    s1.clear()

def test_read_sif():
    s = SIF()
    s.read_sif(simple)
    try:
        s.read_sif("rubbusish")
        assert False
    except:
        assert True

def test_save():
    s1 = SIF(simple)
    import tempfile
    f = tempfile.NamedTemporaryFile(suffix=".sif")
    s1.save(f.name)
    s3 = SIF(f.name)
    assert s1 == s3
    f.close()

def test_sifreader_wrong():
    try:
        sif = SIF("wrong.sif")
        assert False
    except:
        assert True

    try:
        sif = SIF("missing.sif")
        assert False
    except IOError:
        assert True

def test_to_reactions():
    from cno import cnodata
    filename = cnodata("PKN-ToyMMB.sif")
    s = SIF(filename)
    s.add_reaction("a=and1")
    s.add_reaction("b=and1")
    s.add_reaction("and1=c")
    r = s.to_reactions()
    
def test_operators():
    s1 = SIF()
    s1.add_reaction("a=b")

    s2 = SIF()
    s2.add_reaction("a=b")

    s3 = SIF()
    s3.add_reaction("a=c")

    s4 = SIF()
    s4.add_reaction("!a=b")
    assert s1==s2
    assert (s1 == 1) == False

    s2.add_reaction("a=c")
    assert (s1 == s2) == False

    assert (s1 == s4) == False
    
    

    try:
        s2.add_reaction("a")
        assert False
    except:
        assert True

def test_exportsbml():
    s1 = SIF()
    s1.add_reaction("A=C")
    s1.add_reaction("!B=C")
    s1.to_sbmlqual("test.xml")

    s2 = SIF()
    s2.read_sbmlqual("test.xml")
    assert s1 == s2

def test_constructor():
    # bad format
    import tempfile
    fh = tempfile.NamedTemporaryFile(delete=False)
    fh.write("A 1B\n")
    fh.close()
    try:
        s = SIF(fh.name)
        assert False
    except:
        assert True
    fh.delete = True
    fh.close()

    # bad instance
    try:
        s = SIF([1,1,2])
        assert False
    except:
        assert True

def test_remove_and_gates():
    s = SIF()
    s.add_reaction("A+B=1")
    s.add_reaction("A^B=1")
    s.remove_and_gates()

def _test_constructor2():
    """test the convert_ands parameter"""
    fh = tempfile.NamedTemporaryFile(delete=False)
    fh.write("A^B 1 C\n")
    fh.write("a 1 and1\n")
    fh.write("b 1 and1\n")
    fh.write("and1 1 c\n")
    fh.close()

    s = SIF(fh.name, convert_ands=False)
    s = SIF(fh.name, convert_ands=True)

    fh.delete = True
    fh.close()

def test_constructor3():
    """edge can be only 1 or -1 if fomrat is cno  """
    fh = tempfile.NamedTemporaryFile(delete=False)
    fh.write("A activate B\n")
    fh.close()
    try:
        s = SIF(fh.name)
        assert False
    except:
        assert True
    s = SIF(fh.name, frmt='generic')
    fh.delete = True
    fh.close()

def test_constructor_bad_and_gate():
    fh = tempfile.NamedTemporaryFile(delete=False)
    fh.write("A 1 and1\n")
    fh.write("B 1 and1\n")
    fh.write("and1 -1 C\n") # this is not possible
    fh.close()

    try:
        s = SIF(fh.name)
        assert False
    except:
        assert True

    fh.delete = True
    fh.close()

def _test_constructor_ignore_and():
    fh = tempfile.NamedTemporaryFile(delete=False)
    fh.write("A 1 and1\n")
    fh.write("B 1 and1\n")
    fh.write("and1 1 C\n")
    fh.write("a 1 b\n")
    fh.close()

    s = SIF(fh.name, ignore_and=True)

    fh.delete = True
    fh.close()

def test_equal():

    s1 = SIF()
    s1.add_reaction("A=C")

    s2 = SIF()
    s2.add_reaction("A=B")
    assert s1 != s2

    s2 = SIF()
    s2.add_reaction("!A=C")
    assert s1 != s2

def test_plot():
    s = SIF()
    s.add_reaction("A=B")
    s.plot()


def test_notedge():
    s = SIF()
    s.add_reaction("a=b")
    assert s.notedge("-1") == "!"
    assert s.notedge("1") == ""

def test_to_sif():
    fh = easydev.TempFile()
    s1 = SIF(getdata('PKN-test_all.sif'))
    s1.save(fh.name)

    s2 = SIF(fh.name)
    fh.delete()

    assert s1 == s2




