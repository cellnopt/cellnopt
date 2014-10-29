"""This test file cover the cno.midas module"""
from cno.io import Measurement, Measurements, MIDASBuilder
import numpy

from easydev.easytest import assert_list_almost_equal, TempFile

filenames = ['MD-ToyMMB_bis.csv',
             'MD-LiverDREAM.csv']


def test_experiment():
    e = Measurement("AKT", 0, {'TGFa':1}, {'AKT':1}, 10.1)
    e.units = "second"
    try:
        e.units = "dummy"
        assert False
    except:
        assert True
    e.cellLine = "TEST"
    e.protein_name = "zap70"

    e.stimuli = {"TGFa":0}
    assert e.stimuli == {"TGFa":0}

    e.get_cues()
    e.cues_as_dict()

    print(e)



def test_midasbuilder():
    m = MIDASBuilder()
    e1 = Measurement("AKT", 0, {"EGFR":1}, {"AKT":0}, 0.1)
    e2 = Measurement("AKT", 5, {"EGFR":1}, {"AKT":0}, 0.5)
    e3 = Measurement("AKT", 10, {"EGFR":1}, {"AKT":0}, 0.9)
    e4 = Measurement("AKT", 0, {"EGFR":0}, {"AKT":0}, 0.1)
    e5 = Measurement("AKT", 5, {"EGFR":0}, {"AKT":0}, 0.1)
    e6 = Measurement("AKT", 10, {"EGFR":0}, {"AKT":0}, 0.1)
    for e in [e1,e2,e3,e4,e5,e6]:
        m.add_measurement(e)
    #m.to_midas("test.csv")
    




def test_experimentS():

     es = Measurements()
     e1 = Measurement("AKT", 0, {"EGFR":1}, {"AKT":0}, 0.1)
     e2 = Measurement("AKT", 5, {"EGFR":1}, {"AKT":0}, 0.5)
     es.add_single_measurements([e1,e2])
     assert len(es)==2













