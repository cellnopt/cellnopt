from cno.boolean import Steady
from cno import cnodata
from numpy import array
import numpy
from nose.tools import assert_almost_equal

def test_steady():


    # simulation using all bit strings to 1
    s = Steady(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
    s.simulate() 
    assert_almost_equal(s.score(), 0.30906031746, places=10)


    s = Steady(cnodata("PKN-LiverDREAM.sif"), cnodata("MD-LiverDREAM.csv"))
    s._params['include_time_zero'] = False
    s.simulate() 
    assert_almost_equal(s.score(), 0.2574947901818182, places=10)
            
    s._params['include_time_zero'] = True
    s.simulate() 
    assert_almost_equal(s.score(), 0.4696160023030303, places=10)

    # from CNOR. checked 23/4/2015
    s = Steady(cnodata("PKN-LiverDREAM.sif"), cnodata("MD-LiverDREAM.csv"))
    s.simulate() 
    s._params['NAFac'] = 0
    assert_almost_equal(s.score(), 0.4696160023030303, places=10)

    s = Steady(cnodata("PKN-ExtLiverPCB.sif"), cnodata("MD-ExtLiverPCB.csv"))
    s._params['include_time_zero'] = False
    s.simulate() 
    s._params['NAFac'] = 0
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


    # checked against cellnoptr 1.11.5
    s = Steady(cnodata("PKN-LiverDREAM.sif"),
            cnodata("MD-LiverDREAM.csv"))
    s.preprocessing(compression=False)
    s.simulate()
    s.debug_score = True
    assert_almost_equal(s.score(), 0.4696160023030303, places=10)



def test_optimise():
    s = Steady(cnodata('PKN-ToyMMB.sif'), cnodata('MD-ToyMMB.csv'))
    #s.preprocessing()
    s.optimise(verbose=False)
    assert_almost_equal(s.results.results.best_score[-1], 0.0296702380952, places=10)


def test_shift_simulation():
    # In CellNOptR and cno, the simulation stops if the difference between
    # current sim and previous are small or the number of tic is larger than 
    # nSpecies*1.2
    s = Steady(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
    s._shift = -1
    s.simulate()
    assert_almost_equal(s.score(), 0.28835194997119434, places=10)





