from cno.io.obs import OBS
from cno import cnodata, XMIDAS, CNOError
from easydev import TempFile
import numpy as np

def test_obs():
    # Create some data
    midas = XMIDAS(cnodata("MD-ToyMMB.csv"))

    midas.df.loc[('mock', 'experiment_0', 10),'Akt'] = -1

    midas.df.loc[('mock', 'experiment_1', 10),'Akt'] = np.nan

    # testing constructor and accessors
    o = OBS(midas)
    o.midas
    o.midas = midas 

    o = OBS()
    # nothing happends here because obs has no midas yet
    o.save('test.obs')

    try:
        o.set_time(10)
        assert False
    except CNOError:
        assert True

    o.midas = midas

    o.set_time()#fetches automatically the time
    # bad time
    try:
        o.set_time(100)
        assert False
    except CNOError:
        assert True
    # specify the time
    o.set_time(10)

    fh = TempFile()
    o.save(fh.name)
    fh.delete()
