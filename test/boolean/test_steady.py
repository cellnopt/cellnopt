from cno.boolean import Steady
from cno import cnodata
from numpy import array
import numpy

def test_steady():


    # simulation using all bit strings to 1
    s = Steady(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
    s.simulate(1) 
    assert s.score() == 0.094774603174603178


    s = Steady(cnodata("PKN-LiverDREAM.sif"), cnodata("MD-LiverDREAM.csv"))
    s.simulate(1) 
    assert s.score() == 0.2574947901818182

    s = Steady(cnodata("PKN-ExtLiverPCB.sif"), cnodata("MD-ExtLiverPCB.csv"))
    s.simulate(1) 
    assert s.score() == 0.29198998104079721


    # Now test if several reactions are removed
    s = Steady(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
    s.simulate(reactions=['EGF=Ras'])
    assert s.score() == 0.16610942460317463
    s.simulate(reactions=['EGF=Ras', 'EGF=PI3K'])
    assert s.score() == 0.16611567460317461


