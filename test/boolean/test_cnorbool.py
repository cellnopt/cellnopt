from cno.boolean.cnorbool import CNORbool
from cno import cnodata


from nose.plugins.attrib import attr

@attr('skip_travis')
def test_cnorbool():
    c = CNORbool(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
    c.optimise(maxgens=3, popsize=10)
    c.create_report()



