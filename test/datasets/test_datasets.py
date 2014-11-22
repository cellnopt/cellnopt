import os
from cno.datasets import ToyMMB, ToyPB, ToyPCB, ExtLiverPCB
from cno.datasets import cnodata, register

from nose.plugins.attrib import attr


def test_registered_cnodata():
    for filename in register.registered:
        fullpath = cnodata(filename)
        assert os.path.exists(fullpath), (filename,fullpath)

    try:
        cnodata('dummy')
        assert False
    except:
        assert True

def test_datasets_toymmb():
    ToyMMB.model
    ToyMMB.data
    ToyMMB.plot()

@attr('skip')
def test_all_plots():
    import cno.datasets as ds

    ds.ToyPB.plot()
    ds.ToyPCB.plot()
    ds.ExtLiverPCB.plot()
    ds.LiverDREAM.plot()
    ds.ToyPB_True.plot()

def test_cnodata():
    cnodata()
    cnodata('*SBML*')
    cnodata('unknown')
    
