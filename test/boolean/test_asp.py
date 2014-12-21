from cno.boolean import CASPO
from cno import cnodata
from nose.plugins.attrib import attr


@attr('skip_travis')
class TestCASPO(object):
    @classmethod
    def setup_class(klass):
        klass.pkn = cnodata("PKN-ToyMMB.sif")
        klass.midas = cnodata("MD-ToyMMB.csv")
        klass.caspo = CASPO(verbose=True)
        klass.caspo.optimise(klass.pkn, klass.midas, fit=0.5, size=10)


    def _test_plot(self):
        self.caspo.hist_mse()
        self.caspo.hist_model_size()


