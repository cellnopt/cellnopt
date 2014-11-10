"""This test file cover the cno.midas module"""
from cno.io import Measurement, Measurements, MIDASBuilder
import numpy

from easydev.easytest import assert_list_almost_equal, TempFile

filenames = ['MD-ToyMMB_bis.csv',
             'MD-LiverDREAM.csv']


def test_experiment():
    e = Measurement("AKT", 0, {'TGFa':1}, {'AKT':1}, 10.1)
    e.units = "second"
    e.units
    try:
        e.units = "dummy"
        assert False
    except:
        assert True
    e.cellLine = "TEST"
    e.protein_name = "zap70"
    e.time = 10

    e.stimuli = {"TGFa":0}
    assert e.stimuli == {"TGFa":0}
    try:
        e.stimuli = {"TGFa":10}
        assert False
    except:
        assert True

    e.inhibitors = {"AKT":0}
    try:
        e.inhibitors = {"AKT":10}
        assert False
    except:
        assert True
    e.data = 10.

    e.get_cues()
    e.cues_as_dict()

    print(e)


def test_midasbuilder():
    m = MIDASBuilder()

    # what would happen if we try to get XMIDAS from 0 measuremnts
    xm = m.xmidas
    assert len(xm.df) == 0


    e1 = Measurement("AKT", 0, {"EGFR":1}, {"AKT":0}, 0.1)
    e2 = Measurement("AKT", 5, {"EGFR":1}, {"AKT":0}, 0.5)
    e3 = Measurement("AKT", 10, {"EGFR":1}, {"AKT":0}, 0.9)
    e4 = Measurement("AKT", 0, {"EGFR":0}, {"AKT":0}, 0.1)
    e5 = Measurement("AKT", 5, {"EGFR":0}, {"AKT":0}, 0.1)
    e6 = Measurement("AKT", 10, {"EGFR":0}, {"AKT":0}, 0.1)
    m.add_measurements([e1, e2, e3, e4, e5, e6])
    m.test_example()
    assert len(m)>10
    m.get_colnames()
    m.stimuli
    m.inhibitors
    xm = m.xmidas

    # Constructor from existing measurements:



def test_experimentS():

     es = Measurements()
     e1 = Measurement("AKT", 0, {"EGFR":1}, {"AKT":0}, 0.1)
     e2 = Measurement("AKT", 5, {"EGFR":1}, {"AKT":0}, 0.5)
     es.add_measurements([e1,e2])
     assert len(es)==2
     assert es.species == ['AKT']


     m = Measurements(es)
     #assert m.xmidas == xm









