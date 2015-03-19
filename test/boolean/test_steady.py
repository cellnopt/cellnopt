from cno.boolean import Steady
from cno import cnodata
from numpy import array
import numpy
from nose.tools import assert_almost_equal

def test_steady():


    # simulation using all bit strings to 1
    s = Steady(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
    s.simulate(1) 
    assert_almost_equal(s.score(), 0.0947746031746031, places=10)


    s = Steady(cnodata("PKN-LiverDREAM.sif"), cnodata("MD-LiverDREAM.csv"))
    s.simulate(1) 
    assert_almost_equal(s.score(), 0.257494790181818, places=10)

    s = Steady(cnodata("PKN-ExtLiverPCB.sif"), cnodata("MD-ExtLiverPCB.csv"))
    s.simulate(1) 
    assert_almost_equal( s.score(), 0.2919899810407972, places=10)


    # Now test if several reactions are removed
    s = Steady(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
    s.simulate(reactions=['EGF=Ras'])
    assert_almost_equal(s.score(), 0.1661094246031746,places=10)

    s.simulate(reactions=['EGF=Ras', 'EGF=PI3K'])
    assert_almost_equal(s.score(), 0.1661156746031746, places=10)

def test_steady_with_and_gates():
    s = Steady(cnodata('PKN-ToyMMB.sif'), cnodata('MD-ToyMMB.csv'))
    s.preprocessing()
    reactions = ['EGF=PI3K', 'EGF=Raf', 'Mek=Erk', 'Mek=p90RSK', 
            'PI3K=Akt', 'Raf=Mek', 'TNFa=Hsp27',  'TNFa=NFkB', 'TNFa=PI3K']

    s.simulate(reactions=reactions)
    assert_almost_equal(s.score(), 0.02964260651629073, places=10)

