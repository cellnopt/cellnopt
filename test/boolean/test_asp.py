from cno.boolean import CASPO
from cno import cnodata



class TestCASPO(object):
    @classmethod
    def setup_class(klass):
        klass.pkn = cnodata("PKN-ToyMMB.sif")
        klass.midas = cnodata("MD-ToyMMB.csv")
        klass.caspo = CASPO(verbose=True)
        klass.caspo.optimise(klass.pkn, klass.midas, fit=0.5, size=10)


    def test_plot(self):
        self.caspo.hist_mse()
        self.caspo.hist_model_size()


t = TestCASPO()
t.setup_class()
t.test_plot()
