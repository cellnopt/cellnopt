# -*- python -*-
#
#  This file is part of CNO software
#
#  Copyright (c) 2013-2014 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: http://github.com/cellnopt/cellnopt
#
##############################################################################
import sys
import os

import pandas as pd
import pylab

from cno.core import CNOBase, CNORBase
from cno.core.params import  OptionsBase

from cno.core.results import ODEResults
from cno.core import ReportODE
from cno.core.params import ParamsSSM

from biokit.rtools import bool2R

from cno.core.params import params_to_update


__all__ = ["CNORode"]


class CNORode(CNOBase, CNORBase):
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
    def __init__(self, model=None, data=None, tag=None, verbose=True,
                 verboseR=False, config=None, use_cnodata=False):
        CNOBase.__init__(self,model, data, verbose=verbose, tag=tag, 
                config=config, use_cnodata=use_cnodata)
        CNORBase.__init__(self, verboseR=verboseR)

        self._report = ReportODE()
        self._report.Rdependencies = [] ## just to speed up report

        self.results = ODEResults()

        self.config.General.pknmodel.value = self.pknmodel.filename
        self.config.General.data.value = self.data.filename

        p = ParamsSSM()
        self.config.add_section(p)

        self._library = 'CNORode'
        #CNORodePBMstNeu

    def _init(self):
        script_template = """
        library(%(library)s)
        pknmodel = readSIF("%(pknmodel)s")
        cnolist = CNOlist("%(midas)s")
        model = preprocessing(cnolist, pknmodel, compression=%(compression)s,
                   expansion=%(expansion)s, cutNONC=%(cutnonc)s,
                   maxInputsPerGate=%(maxInputsPerGate)s)
        reactions = model$reacID
        species = colnames(cnolist@signals[[1]])"""
        params = {
                'library': self._library,
                'pknmodel': self.pknmodel.filename,
                'midas': self.data.filename,
            'compression': bool2R(self._compression),
            'expansion': bool2R(self._expansion),
            'cutnonc': bool2R(self._cutnonc),
             'maxInputsPerGate': self._max_inputs_per_gate,
                }
        self.session.run(script_template % params)
        self.species = self.session.species

    @params_to_update
    def optimise(self,  n_diverse=10, dim_ref_set=10, maxtime=60,
                 verbose=False, reltol=1e-4, atol=1e-3, maxeval='Inf',
                 transfer_function=3, maxstepsize='Inf', reuse_ode_params=False,
                 local_solver=None):
        """Optimise the ODE parameters using SSM algorithm 

        :param int maxtime: (default 10)
        :param int ndiverse: (default 10)
        :param int dim_refset: (default 10)
        :param bool: ode_params if True, load the ode_params.RDAta file saved in a previous run
            mus be compatible with the model.
        
        verbose should be False all the time internally to the R code. Here, verbose
        meana we want to see the status of the optimisation (not all warnings and R
        errors).
        """
        self.logging.info("Running the optimisation. Can take a very long"
                          "time. To see the progression, set verboseR "
                          "attribute to True")
        # update config GA section with user parameters
        self._update_config('SSM', self.optimise.actual_kwargs)
        ssmd = self.config.SSM.as_dict()

        if self.session.get('ode_params') is None:
            self.session.run('ode_params=NULL')
        if reuse_ode_params is False:
            self.session.run('ode_params=NULL')

        # todo: ode_params to be provided as input
        script = """
        library(%(library)s)
        pknmodel = readSIF("%(pknmodel)s")
        cnolist = CNOlist("%(midas)s")


        model = preprocessing(cnolist, pknmodel, compression=%(compression)s,
                   expansion=%(expansion)s, cutNONC=%(cutnonc)s,
                   maxInputsPerGate=%(maxInputsPerGate)s)


        reactions = model$reacID
        species = colnames(cnolist@signals[[1]])


        if (is.null(ode_params) == TRUE){
         ode_params = createLBodeContPars(model)
        }
        ode_params = parEstimationLBodeSSm(cnolist, model, 
            maxtime=%(maxtime)s, maxStepSize=%(maxstepsize)s, dim_refset=%(dim_ref_set)s, maxeval=%(maxeval)s,
            verbose=F, ndiverse=%(n_diverse)s, ode_parameters=ode_params,
            local_solver=%(local_solver)s)
        """

        if local_solver is None:
            local_solver = 'NULL'
        params = {
                'library': self._library,
            'pknmodel': self.pknmodel.filename,
            'midas': self.data.filename,
            'compression': bool2R(self._compression),
            'expansion': bool2R(self._expansion),
            'cutnonc': bool2R(self._cutnonc),
             'maxInputsPerGate': self._max_inputs_per_gate,
             'local_solver': local_solver
            }

        params.update(ssmd)

        self.session.run(script % params)

        ssm_results = self.session.ode_params['ssm_results'].copy()
        self.ssm_results = ssm_results
        results = {
                'best_score': ssm_results['fbest'],
                'all_scores': ssm_results['f'],
                'reactions': self.session.reactions[:],
                'best_params': ssm_results['xbest']
        }
        for k,v in self.session.ode_params.items():
            results[k] = v.copy()

        self.results = ODEResults()
        self.results.results = results
        self.species = self.session.species

    def plot_errors(self, show=True):
        self._set_simulation()
        self.midas.plot(mode="mse") 
        #self.midas.plotSim()
        if show is False:
            pylab.close()

    def simulate(self, params, verboseR=False):
        # The first call is slow but then, it is faster but still
        # 10 times slower than the pure R version
        save_verboseR = self.verboseR
        self.verboseR = verboseR
        if self.session.get("simulator_initialised") is None:
            script = """
                library(%(library)s)
                pknmodel = readSIF("%(pknmodel)s")
                cnolist = CNOlist("%(midas)s")
                model = preprocessing(cnolist, pknmodel, compression=%(compression)s,
                   expansion=%(expansion)s, cutNONC=%(cutnonc)s,
                   maxInputsPerGate=%(maxInputsPerGate)s)
                indices = indexFinder(cnolist, model,verbose=FALSE)
                ode_params = createLBodeContPars(model)
                objective_function = getLBodeContObjFunction(cnolist, model,
                    ode_params, indices)
                simulator_initialised = T
            """
            pars = {
                'library': self._library,
                'pknmodel': self.pknmodel.filename,
                'midas': self.data.filename,
            'compression': bool2R(self._compression),
            'expansion': bool2R(self._expansion),
            'cutnonc': bool2R(self._cutnonc),
             'maxInputsPerGate': self._max_inputs_per_gate,
            }
            self.session.run(script % pars)

        self.session['params'] = params
        script = """
            score = objective_function(params)
        """
        self.session.run(script)
        self.verboseR = save_verboseR
        return self.session.score

    def get_sim_data(self, bs=None):
        """

        input could be a bitstring with correct length and same order
        OR a model

        """
        if bs is None:
            bs = self.results.results.parValues
        else:
            # TODO check assert length bs is correct
            pass

        script_template = """
        library(%(library)s)
        pknmodel = readSIF("%(pknmodel)s")


        cnolist = CNOlist("%(midas)s")
        model = preprocessing(cnolist, pknmodel, compression=%(compression)s,
                   expansion=%(expansion)s, cutNONC=%(cutnonc)s,
                   maxInputsPerGate=%(maxInputsPerGate)s)
        sim_data = plotLBodeFitness(cnolist,model, ode_parameters=ode_params)
        """

        params = {
                'library': self._library,
                'pknmodel': self.pknmodel.filename,
                'midas': self.data.filename,

            'compression': bool2R(self._compression),
            'expansion': bool2R(self._expansion),
            'cutnonc': bool2R(self._cutnonc),
             'maxInputsPerGate': self._max_inputs_per_gate,
                #'tag':tag
                }
        script = script_template % params
        self.session.run(script)
        # FIXME what about species/experiments

        sim_data = self.session.sim_data
        self.sim = pd.concat([pd.DataFrame(x, columns=self.species) 
            for x in sim_data])

    def _get_models(self):
        return self.results.cnorbool.models
    models = property(_get_models)

    def _set_simulation(self):
        self.get_sim_data()
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
        r = ODEParameters(self.results.results.parNames, self.results.results.parValues)
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
        r = ODEParameters(self.results.results.parNames, self.results.results.parValues)
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

    def create_report_images(self):
        model = self.cnograph.copy()
        model.plot(filename=self._report._make_filename("pknmodel.svg"), show=False)
        model.preprocessing()
        model.plot(filename=self._report._make_filename("expmodel.png"), show=False)
        self._plot_ode_parameters_k(filename=self._report._make_filename("ode_parameters_k.png"), 
                show=False)
        self._plot_ode_parameters_n(filename=self._report._make_filename("ode_parameters_n.png"), 
                show=False)

        self.plot_errors(show=True)
        self._report.savefig("Errors.png")
        
        self.midas.plot()
        self._report.savefig("midas.png")

        pylab.close()
        self.plot_fitness(show=True, save=False)
        self._report.savefig("fitness.png")

    def plot_fitness(self, show=True, save=False):
        self.results.plot_fit()

        if save is True:
            self._report.savefig("fitness.png")

        if show is False:
            pylab.close()

    def create_report(self):
        self._create_report_header()

        txt = """<pre class="literal-block">\n"""
        #txt += "\n".join([x for x in self._script_optim.split("\n") if "write.csv" not in x])
        txt += "o.report()\n</pre>\n"
        self._report.add_section(txt, "Script used")
        
        txt = """<a href="http://www.cellnopt.org/">
            <object data="pknmodel.svg" type="image/svg+xml">
            <span>Your browser doesn't support SVG images</span> </object></a>"""
        txt += """<a class="reference external image-reference" href="scripts/exercice_3.py">
<img alt="MIDAS" class="align-right" src="midas.png" /></a>"""
        
        self._report.add_section(txt, "PKN graph", [("http://www.cellnopt.org", "cnograph")])
                                                                                            
        self._report.add_section('<img src="expmodel.png">', "Expanded before optimisation")
        self._report.add_section( """
        <img src="expmodel.png">'
        <img src="ode_parameters_k.png">
        <img src="ode_parameters_n.png">
        """, "Optimised model")
               
        self._report.add_section('<img src="fitness.png">', "Fitness")
        self._report.add_section('<img src="Errors.png">', "Errors")

        self._report.add_section(self._report.get_html_reproduce(), "Reproducibility")
        fh = open(self._report.report_directory + os.sep + "rerun.py", 'w')
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
        self._report.add_section(txt, "stats")

        # dependencies
        #table = self.get_table_dependencies()
        #fh.write(table.to_html())

        self._report.write("index.html")

    def _get_stats(self):
        res = {}
        #res['Computation time'] = self.total_time
        try:
            res['Best Score'] = self.results.results['best_score']
        except:
            pass
        return res


class ODEParameters(object):
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


def standalone(args=None):
    """This function is used by the standalone application called cellnopt_boolean

    ::

        cellnopt_ode --help

    """
    if args is None:
        args = sys.argv[:]

    from cno.core.standalone import Standalone
    user_options = OptionsODE()
    stander = Standalone(args, user_options)

    # just an alias
    options = stander.options

    if options.onweb is True or options.report is True:
        trainer = CNORode(options.pknmodel, options.data, verbose=options.verbose,
            verboseR=options.verboseR, config=options.config_file, 
             use_cnodata=options.cnodata)
        trainer.preprocessing()
    else:
        stander.help()

    trainer.optimise(**stander.user_options.config.SSM.as_dict())

    stander.trainer = trainer
    stander.report()


class OptionsODE(OptionsBase):
    def __init__(self):
        prog = "cno_ode_steady"
        version = prog + " v1.0 (Thomas Cokelaer @2014)"
        super(OptionsODE, self).__init__(version=version, prog=prog)
        self.add_section(ParamsSSM())


if __name__ == "__main__":
    """Used by setup.py as an entry point to :func:`standalone`"""
    standalone(sys.argv)










