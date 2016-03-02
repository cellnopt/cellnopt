from cno import cnodata
from cno.feeder import Feeder
from nose.plugins.attrib import attr


@attr('skip_travis')
def test_feeder():
    f = Feeder()
    pknmodel = cnodata("PKN-ToyMMB.sif")
    midas = cnodata("MD-ToyMMB.csv")
    f.run(model=pknmodel, data=midas)
    f.newlinks
    print(f)
    f.plot()

