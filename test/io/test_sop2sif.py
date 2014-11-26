import tempfile
from cno.io.sop2sif import SOP2SIF
from cno import getdata

filename = getdata("test_data.sop") 


def test_sop2sif():
    s2s = SOP2SIF(filename)
    f = tempfile.NamedTemporaryFile()
    s2s.export2sif(filename=f.name, include_and_gates=False)
    f = tempfile.NamedTemporaryFile()
    s2s.export2sif(filename=f.name, include_and_gates=True)

    #assert len(s2s.reacID) == 173
    assert len(s2s.reactions) == 348
    
    assert len([x for x in s2s.species if x.startswith("and")==False])
    len(s2s)



