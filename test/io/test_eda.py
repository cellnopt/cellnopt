from cno import EDA
from cno import SIF
from cno.testing import getdata
import tempfile

filename = getdata("test_simple.eda") 
sifname = getdata("test_simple.sif") 


def test_eda():
    e = EDA(filename, verbose=True)
    s1 = e.export2sif()
    s2 = SIF(sifname)
    assert s1 == s2

    try:
        e = EDA("dummy")
        assert False
    except:
        assert True


def test_eda_format():
    # test eda with header (used in HPN-Dream8 challenge
    fh = tempfile.NamedTemporaryFile(delete=False)
    fh.write("EdgeScore\n")
    fh.write("A (1) B = 0.5\n")
    fh.close()
    e = EDA(fh.name)
    fh.delete = True
    fh.close()

    # test bad data sets

    fh = tempfile.NamedTemporaryFile(delete=False)
    fh.write("A (1) B = 0.5\n")
    fh.write("\n") # should be skipped without errors
    fh.write("A B = 0.5\n") # should raise an erro
    fh.close()
    try:
        e = EDA(fh.name)
        assert False
    except:
        assert True
    fh.delete = True
    fh.close()

