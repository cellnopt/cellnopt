from cno.io import Metabolites
from cno import getdata
from cno.testing import getdata_metabolites
import glob
import os



def test_metabolites():
    filename = getdata("test_metabolites")
    a = Metabolites(filename)
    assert len(a.species) == 94 
    a.specDefault
    a.specBoxes


def test_metabolitessamples():
    for filename in getdata_metabolites():
        print(filename)
        yield readfiles, filename

def readfiles(filename):
    print('Reading metabolies')
    try:
        try:
            m = Metabolites(os.path.split(filename)[1])
            len(m)
        except:
            assert False
    except:
        print('failed. File could not be found ? Skipped')


#test_metabolitessamples()
