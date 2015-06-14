from cno.boolean.cnordt import CNORdt
from cno import cnodata


def test_cnordt():
    c = CNORdt(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
    c.optimise(maxgens=3, popsize=10, bool_updates=2)
    c.create_report()
    # cleanup
    c.cleanup()
