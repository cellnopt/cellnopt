from cno.testing import getdata
from cno import SIF, CNOError

def test_getdata():
    filename = getdata("test_simple.sif")
    s = SIF(filename)

    try:
        filename = getdata("unknown")
        assert False
    except CNOError:
        assert True
    except:
        # should not enter here since a CNOError is raised.
        assert False

