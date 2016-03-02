from cno import cnodata, CNORode

from nose.plugins.attrib import attr


@attr('skip_travis')
def test_cnorode():
    c = CNORode(cnodata("PKN-ToyMMB_T2.sif"), cnodata("MD-ToyMMB_T2.csv"), verboseR=True)
    c.config.SSM.n_diverse.value = 10
    c.config.SSM.maxtime.value = 5
    c.optimise()
    c.plot_ode_parameters()
    c._plot_ode_parameters_n()
    c._plot_ode_parameters_k()
    c._plot_ode_parameters_n()
    c._plot_ode_parameters_k()
    c._plot_ode_parameters_n()
    c._plot_ode_parameters_k()

