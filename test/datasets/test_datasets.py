import os
from cno.datasets import ToyMMB
from cno.datasets import cnodata, registers


def test_registered_cnodata():
    for filename in registers:
        fullpath = cnodata(filename)
        assert os.path.exists(fullpath)

    try:
        cnodata('dummy')
        assert False
    except:
        assert True

def test_datasets_toymmb():

    ToyMMB.model_filename
    ToyMMB.data_filename
    ToyMMB.description
    ToyMMB.plot()
