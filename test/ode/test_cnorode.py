from cno import cnodata, CNORode


def test_cnorode():
    c = CNORode(cnodata("PKN-ToyMMB_T2.sif"), cnodata("MD-ToyMMB_T2.csv"), verboseR=True)
    c.config.SSM.n_diverse = 10
    c.config.SSM.maxtime = 10
    c.optimise()
    c.plot_ode_parameters()
    c._plot_ode_parameters_n()
    c._plot_ode_parameters_k()
    c._plot_ode_parameters_n()
    c._plot_ode_parameters_k()
    c._plot_ode_parameters_n()
    c._plot_ode_parameters_k()

