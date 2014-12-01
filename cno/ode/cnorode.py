import sys
import os

import pandas as pd
import pylab

from easydev import Logging, AttrDict

from cno.io.multigraph import CNOGraphMultiEdges
from cno import CNOGraph, XMIDAS
from cno.core import CNOBase

#from htmltools import HTMLReportODE
#from cnobase import CNObase, OptionBase

from cno.core.params import Parameters
from cno.misc.results import ODEResults

from biokit.rtools.session import RSession
from biokit.rtools import bool2R


__all__ = ["CNORode"]



class ODEParameters(Parameters):
    # THe keys used here have the same caps as in the R code.
    ssm_params = {
        "maxtime": 60,
        "ndiverse": 10,
        "dim_refset": 10,
        "transfer_function": 3
    }
    
    def init_ssm(self):
        default = self.ssm_params   # just an alias
        self.add_parameter("SSM", "--maxtime", "maxtime", default['maxtime'], 
                           "The elitism number (should be 10% of the popsize)")
        self.add_parameter("SSM", "--dim-ref-set", "dim_refset", default['dim_refset'], 
                           "The elitism number (should be 10% of the popsize)")
        self.add_parameter("SSM", "--n-diverse", "ndiverse", default['ndiverse'], 
                           "The elitism number (should be 10% of the popsize)")
        self.add_parameter("SSM", "--ssm-verbose", "verbose", True, 
                           "The elitism number (should be 10% of the popsize)")
   
   



class SimulateODE(object):
    """DRAFT. examlpe of script to generate the sim and save in a file"""     
    def __init__(self, pknmode, midas):
        script = """
        library(CellNOptR)
        pknmodel = readSIF("%(pknmodel)s")
        cnolist = CNOlist("%(midas)s")
        
        load("params.RData")
        sim = plotLBodeFitness(cnolist, pknmodel, res, 
                               initial_state=0.1)
        signals = colnames(cnolist@signals$`0`)
        colnames(output$mse) = signals
        for (i in seq_along(sim[[1]])){
            colnames(sim[[1]][[i]]) = signals
        }

        
        write.csv(sim, "sim.csv")""" % {'pknmodel':pknmodel, 'midas':midas}
        self.runRscript(script)






class CNORode(CNOBase):
    """Access to CellNOptR R package to run boolean analysis

    ::

        c = pipeline.CNObool("PKN-test.sif", "MD-test.csv")
        c.optimise(compression=True, expansion=True, reltol=.15)


    Results are stored in :attr:`results`. Information stored are various.
    The errors corresponding to the best models can be visualised with :meth:`plot_errors`
    and models within the tolerance are stored in :attr:`models.

    .. plot::
        :include-source:

        from cno import cnodata, CNORbool
        c = CNORode(cnodata("PKN-ToyMMB.sif"),
            cnodata("MD-ToyMMB.csv"))
        c.optimise()
        c.plot_errors()

    """
    def __init__(self, model, data, verbose=True, verboseR=False):
        super(CNORode, self).__init__(model, data, verbose=verbose)

        # this code could be moved to CNOBase
        self._verboseR = verboseR

        self.session = RSession(verbose=self.verboseR)

        self.ssm_params = {
            "maxtime": 60,
            "ndiverse": 10,
            "dim_refset": 10,
            "transfer_function": 3}

    def _get_verboseR(self):
        return self._verboseR
    def _set_verboseR(self, value):
        self._verboseR = value
        self.session.dump_stdout = value
    verboseR = property(_get_verboseR, _set_verboseR)

    def optimise(self, tag="cnorode", 
            expansion=True, compression=True):
        """Optimise the ODE parameters using SSM algorithm 

        :param int maxtime: (default 10)
        :param int ndiverse: (default 10)
        :param int dim_refset: (default 10)
        :param bool: ode_params if True, load the ode_params.RDAta file saved in a previous run
            mus be compatible with th emodel.        
        
        verbose should be False all the time internally to the R code. Here, verbose
        meana we want to see the status of the optimisation (not all warnings and R
        errors).
        """

        script_template = """
        library(CNORode)
        pknmodel = readSIF("%(pknmodel)s")
        cnolist = CNOlist("%(midas)s")
        reactions = pknmodel$reacID
        
        ode_params = createLBodeContPars(pknmodel)
        ode_params = parEstimationLBodeSSm(cnolist, pknmodel, 
            maxtime=%(maxtime)s, dim_refset=%(dim_refset)s, 
            verbose=F, ndiverse=%(ndiverse)s, ode_parameters=ode_params)
        """
        # to be retrieved inside Python code
        params = {
                'pknmodel': self.pknmodel.filename,
                'midas': self.data.filename,
                'tag':tag
                }
        for k,v in self.ssm_params.items():
            params[k] = v

        script = script_template % params
        self.session.run(script)

        # need to change type of some columns, which are all string
        #results = self.session.results

        ssm_results = self.session.ode_params['ssm_results'].copy()
        self._results = {
                'best_score': ssm_results['fbest'],
                'all_scores': ssm_results['f'],
                'reactions': self.session.reactions[:],
                'tag': tag,
        }
        for k,v in self.session.ode_params.items():
            self._results[k] = v.copy()

        self.results = ODEResults()
        self.results.add_results(self._results)

    def plot_errors(self, show=False):
        self._set_simulation()
        self.midas.plot(mode="mse") 
        self.midas.plotSim()
        self.savefig("errors.png")
        if show == False:
            pylab.close()

    def simulate(self, bs=None, compression=True, expansion=True):
        """

        input could be a bitstring with correct length and same order
        OR a model

        """
        if bs == None:
            bs = self.results.results.parValues
        else:
            # TODO echk assert length bs is correct
            pass
        script_template = """
        library(CNORode)
        pknmodel = readSIF("%(pknmodel)s")
        cnolist = CNOlist("%(midas)s")
        sim_data = plotLBodeFitness(cnolist, pknmodel, ode_parameters=ode_params)
        """
        #model = preprocessing(cnolist, pknmodel, compression=%(compression)s,
        #    expansion=%(expansion)s, maxInputsPerGate=3)
        #mse = simulateODE(cnolist, model, %(bs)s)

        params = {
                'pknmodel': self.pknmodel.filename,
                'midas': self.data.filename,
                #'tag':tag
                }
        script = script_template % params
        self.session.run(script)
        # FIXME what about species/experiments
        self.sim = pd.dataFrame(self.session.sim_data)
        return self.session.sim_data

    def _get_models(self):
        return self.results.cnorbool.models
    models = property(_get_models)


    def _set_simulation(self):
        self.simulate()
        self.midas.create_random_simulation()

        Ntimes = len(self.midas.times)
        Nexp = len(self.midas.experiments.index)
        sim = self.sim.copy()
        
        sim['time'] = [time for time in self.midas.times for x in range(0, Nexp)]
        sim['experiment'] = list(self.midas.experiments.index) * Ntimes
        sim['cellLine'] = [self.midas.cellLines[0]] * sim.shape[0]
        sim.set_index(['cellLine', 'experiment', 'time'], inplace=True)
        sim.sortlevel(1, inplace=True)

        self.midas.sim = sim.copy()


    def plot_ode_parameters(self, **kargs):
        pylab.figure(1);
        self._plot_ode_parameters_k(**kargs)
        pylab.figure(2)
        self._plot_ode_parameters_n(**kargs)

    def _plot_ode_parameters_k(self, **kargs):
        kargs["edge_attribute"] = "ode_k"
        r = ODEParameters2(self.results.results.parNames, self.results.results.parValues)
        data = r.get_k()
        for e in self.cnograph.edges():
            try:
                self.cnograph.edge[e[0]][e[1]]["ode_k"] = data[e[0]][e[1]]
                self.cnograph.edge[e[0]][e[1]]["label"] = " k=%.2f" % data[e[0]][e[1]]
            except:
                self.cnograph.edge[e[0]][e[1]]["ode_k"] = -1
                self.cnograph.edge[e[0]][e[1]]["label"] = " k=??"
                print(e)
        self.cnograph.plot(**kargs)
        # cleanup the label
        for e in self.cnograph.edges():
            del self.cnograph.edge[e[0]][e[1]]["label"]

    def _plot_ode_parameters_n(self, **kargs):
        kargs["edge_attribute"] = "ode_n"
        r = ODEParameters2(self.results.results.parNames, self.results.results.parValues)
        data = r.get_n()
        for e in self.cnograph.edges():
            try:
                self.cnograph.edge[e[0]][e[1]]["ode_n"] = data[e[0]][e[1]]
                self.cnograph.edge[e[0]][e[1]]["label"] = " n=%.2f" % data[e[0]][e[1]]
            except:
                self.cnograph.edge[e[0]][e[1]]["ode_n"] = -1
                self.cnograph.edge[e[0]][e[1]]["label"] = " n=??"
                print(e)
        self.cnograph.plot(**kargs)
        # cleanup the label
        for e in self.cnograph.edges():
            del self.cnograph.edge[e[0]][e[1]]["label"]

    def _update_config(self):
        self.config.add_section("SSM")
        self.config.add_option("SSM", "maxtime", params["maxtime"])
        self.config.add_option("SSM", "n-diverse", params["ndiverse"])
        self.config.add_option("SSM", "dim-ref-set", params["dim_refset"])
        self.config.add_option("SSM", "verbose", True)

    def create_report_images(self):
        self._pknmodel.plot(filename=self._make_filename("pknmodel.png"), show=False)
        self.cnograph.plot(filename=self._make_filename("expmodel.png"), show=False)
        self._plot_ode_parameters_k(filename=self._make_filename("ode_parameters_k.png"), show=False)
        self._plot_ode_parameters_n(filename=self._make_filename("ode_parameters_n.png"), show=False)
        self.plot_errors(show=False)
        self.midas.plot()
        self.savefig("midas.png")
        pylab.close()
        self.plot_fitness(show=False, save=True)

    def report(self, filename="index.html", directory="report", browse=True):
        """Creates the boolean report
        
        Should not create any image here only the report
        """
        self.directory = self._init_report()
        self.create_report_images()
            
        report = HTMLReportODE(self)        
        # complete the report
        self._report(report)
        if browse:
            from browse import browse as bs
            bs(report.directory + os.sep + filename)
        self.save_config_file()

    def _report(self, report):
        directory = self._init_report()
        
        # Save filenames and report in a section
        fname = directory + os.sep + "PKN-pipeline.sif"
        self.cnograph.export2sif(fname)        
        fname = directory + os.sep + "MD-pipeline.csv"
        self.midas.save2midas(fname)

        txt = '<ul><li><a href="PKN-pipeline.sif">input model (PKN)</a></li>'        
        txt += '<li><a href="MD-pipeline.csv">input data (MIDAS)</a></li>'
        txt += '<li><a href="config.ini">Config file</a></li>'
        txt += '<li><a href="rerun.py">Script</a></li></ul>'
        txt += "<bold>some basic stats about the pkn and data e.g. number of species ? or in the pkn section?</bold>"
        
        report.add_section(txt, "Input data files")
        report.add_section(
        """       
         <div class="section" id="Script_used">
         <object height=120 width=300 type='text/x-scriptlet' border=1
         data="description.html"></object>
         </div>""", "Description")

        txt = """<pre class="literal-block">\n"""
        txt += "\n".join([x for x in self._script_optim.split("\n") if "write.csv" not in x])
        txt += "o.report()\n</pre>\n"
        report.add_section(txt, "Script used")
        
        txt = """<a href="http://www.cellnopt.org/">
            <object data="pknmodel.svg" type="image/svg+xml">
            <span>Your browser doesn't support SVG images</span> </object></a>"""
        txt += """<a class="reference external image-reference" href="scripts/exercice_3.py">
<img alt="MIDAS" class="align-right" src="midas.png" /></a>"""
        
        report.add_section(txt, "PKN graph", [("http://www.cellnopt.org", "cnograph")])
                                                                                            
        report.add_section('<img src="expmodel.png">', "Expanded before optimisation")
        report.add_section( """
        <img src="expmodel.png">'
        <img src="ode_parameters_k.png">
        <img src="ode_parameters_n.png">
        """, "Optimised model")
        
        report.add_section('<img src="errors.png">', "Errors")

        report.add_section(self.get_html_reproduce(), "Reproducibility")
        fh = open(self.report_directory + os.sep + "rerun.py", 'w')
        fh.write("from cellnopt.pipeline import *\n")
        fh.write("CNOode(config=config.ini)\n")
        fh.write("c.optimise()\n")
        fh.write("c.report()\n")
        fh.close()

        # some stats
        stats = self._get_stats()
        txt = "<table>"
        for k,v in stats.iteritems():
            txt += "<tr><td>%s</td><td>%s</td></tr>" % (k,v)
        txt += "</table>"
        txt += """<img id="img" onclick='changeImage();' src="fit_over_time.png">\n"""
        report.add_section(txt, "stats")

        # dependencies
        #table = self.get_table_dependencies()
        #fh.write(table.to_html())

        report.write(self.report_directory, "index.html")













class ODEParameters2(object):
    """A class to handle ODE parameters returned by the R package

    """
    def __init__(self, parNames, parValues):
        self.parValues = parValues
        self.parNames = parNames

    def get_n(self):
        res = {}
        for k,v in zip(self.parNames, self.parValues):
            if "_n_" in k:
                e1, e2 = k.split("_n_")
                if e1 not in res.keys():
                    res[e1] = {}
                if e2 not in res[e1].keys():
                    res[e1][e2] = v
                else:
                    raise KeyError
        return res

    def get_k(self):
        res = {}
        for k,v in zip(self.parNames, self.parValues):
            if "_k_" in k:
                e1, e2 = k.split("_k_")
                if e1 not in res.keys():
                    res[e1] = {}
                if e2 not in res[e1].keys():
                    res[e1][e2] = v
                else:
                    raise KeyError
        return res
    def get_tau(self):
        res = {}
        for k,v in zip(self.parNames, self.parValues):
            if "_n_" in k:
                k = k.strip("tau_")
                res[k] = v
        return res

