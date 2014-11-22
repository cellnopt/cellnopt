from cno.io.mapback import MapBack
from cno import CNOGraph, cnodata


def test_mapback():

    # init
    pknmodel = CNOGraph(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
    model = CNOGraph(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
    model.preprocessing()
    mp = MapBack(pknmodel, model)



    links2map = ['Erk^TNFa=Hsp27', '!Akt=Mek']
    new =  mp.mapback(links2map)
    assert sorted(new) == sorted(['!Akt=Mek', 'Akt=Mek', 'Erk=Hsp27', 
        'TNFa=TRAF6', 'TRAF6=p38', 'p38=Hsp27'])



    mp.plotall(links2map)




def test_mapback_and_gates_in_pkn():


    model = CNOGraph(cnodata("PKN-ToyPB_True.sif"), cnodata("MD-ToyPB_True.csv"))
    model.preprocessing(expansion=False)
    pknmodel = CNOGraph(cnodata("PKN-ToyPB_True.sif"), cnodata("MD-ToyPB_True.csv"))
    mp = MapBack(pknmodel, model)

    solutions = sorted(['!erk=ph', '!ikb=nfkb', 'egf=egfr', 'mkk4=p38',
                'ras=map3k1', 'sos=ras', 'tnfa=tnfr', 'tnfr=traf2',
                'traf2=map3k7', 'map3k1^map3k7=mkk4', 'egfr^!ph=sos'])
    new = mp.mapback(['sos^tnfa=p38', '!ikb=nfkb', 'egf^!erk=sos'])
    assert sorted(new) == sorted(solutions)

