from cno import cnodata
from cno.feeder import Feeder


def test_feeder():
    f = Feeder()
    pknmodel = cnodata("PKN-ToyMMB.sif")
    midas = cnodata("MD-ToyMMB.csv")
    f.run(model=pknmodel, data=midas)
    f.newlinks
    print(f)
    f.plot()

