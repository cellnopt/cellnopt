from cno.io.sif import SIF
from cno import cnodata
from easydev import TempFile


def test_sbmlqual():

    # a simple model
    s1 = SIF()
    s2 = SIF(cnodata("PKN-ToyMMB.sif"))
    fh = TempFile()
    s2.to_sbmlqual(fh.name)
    s1.read_sbmlqual(fh.name)
    fh.delete()
    assert s1 == s2



    # a simple model with AND gates
    s1 = SIF()
    s2 = SIF(cnodata("PKN-ToyPB.sif"))
    fh = TempFile()
    s2.to_sbmlqual(fh.name)
    s1.read_sbmlqual(fh.name)
    fh.delete()
    assert s1 == s2



