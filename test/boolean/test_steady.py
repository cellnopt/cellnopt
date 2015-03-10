from cno.boolean import Steady
from cno import cnodata
from numpy import array
import numpy

def test_steady():

    s = Steady(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
    s.simulate(1) 
    assert s.score() == 0.094774603174603178


    # ???
    #s = Steady(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
    #s.preprocessing()
    #reactions = ['!Akt=Mek', '!Akt^Raf=Mek', 'EGF=PI3K', 'Erk=Hsp27',
    #         'Jnk=cJun','Mek=Erk', 'Mek=p90RSK', 'PI3K=Akt',
    #            'Raf=Mek',  'TNFa=Hsp27', 'TNFa=Jnk', 'TNFa=NFkB', 'TNFa=PI3K']
    #s.simulate(reactions)
    #s.score()/2. == 0.10832539682539682


    s = Steady(cnodata("PKN-LiverDREAM.sif"), cnodata("MD-LiverDREAM.csv"))
    s.simulate(1) 
    assert s.score() == 0.2574947901818182





