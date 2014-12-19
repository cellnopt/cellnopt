from cno.core.base import CNORBase, CNOBase
from cno import cnodata
from easydev import TempFile


# To test some of the base functions, need to use something else such as cnorbool

def test_cnobase():

    c = CNOBase(cnodata('PKN-ToyMMB.sif'), cnodata("MD-ToyMMB.csv"))
    c.pknmodel
    c.midas
    c.data
    c.preprocessing()
    c.plot_pknmodel()

    assert c._reac_cnor2cno(['A+B=C']) == ['A^B=C']

    c.plot_midas()
    c.plot_midas(xkcd=True)
    c.config

    fh = TempFile()
    c.save_config_file(fh.name)
    
    c = CNOBase(cnodata('PKN-ToyMMB.sif'), cnodata("MD-ToyMMB.csv"), config=fh.name)

    try:
        c.create_report()
        assert False
    except:
        assert True

    try:
        c.create_report_images()
        assert False
    except:
        assert True


from cno.boolean.cnorbool import CNORbool
def test_cnobase_with_cnorbool():
    c = CNORbool(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"), verbose=True)
    c.verboseR = True
    c.verboseR = False
    c.verbose = False
    c.optimise(maxgens=5, popsize=10) 
    c.plot_fitness()
    c.plot_model()
    c.plot_optimised_model()
    c.plot_mapback_model()
    c._create_report_header()
    c.onweb()











