# -*- python -*-
#
#  This file is part of CNO software
#
#  Copyright (c) 2014 - EBI-EMBL
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
import os
import sys

import numpy as np
import pandas as pd
import pylab

from cno.core import CNOBase
from cno.core.results import BooleanResults
from cno.core.models import BooleanModels
from cno.core import ReportBool
from cno.core.params import OptionsBase, ParamsGA, ParamsGA2


from cno.core.params import params_to_update

import easydev

__all__ = ["CNObool"]


class CNObool(CNOBase):
    """A Python implementation of binary/boolean analysis

    ::

        from cno import CNObool, cnodata
        c = CNObool(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"))
        c.optimise(compression=True, expansion=True, reltol=.15)


    Results are stored in :attr:`results`. Information stored are various.
    The errors corresponding to the best models can be visualised with
    :meth:`plot_errors`  and models within the tolerance are stored in
    :attr:`models`.

    .. plot::
        :include-source:

        from cno import cnodata, CNObool
        c = CNObool(cnodata("PKN-ToyMMB.sif"),
            cnodata("MD-ToyMMB.csv"))
        c.optimise()
        c.plot_errors()

    If you have 2 time points in addition to time zero, the second
    time point can also be optimised::

         c.optimise2()

    Calling plot_errors again will show the 2 time points results.

    .. note:: When calling :meth:`optimise` or :meth:`optimise2` the
        sections of the configuration files called GA and GA2 are updated
        with the parameters provided by the user.
    """
    # params = BooleanParameters.default
    def __init__(self, model, data, tag=None, verbose=True,
                 verboseR=False, config=None, use_cnodata=False):
        """.. rubric:: Constructor

        :param model: model in SIF or SBMLqual format
        :param data: a MIDAS file
        :param tag: not yet used
        :param config: a configuration file stored in :attr:`config`
        :param bool verboseR: switch on/off verbosity of the R session
        :return:

        """
        CNOBase.__init__(self, model, data, tag=tag, verbose=verbose,
                         config=config, use_cnodata=use_cnodata)
        self.logging.info("Initialise R session")
        CNOBase.__init__(self, verboseR)

        self._report = ReportBool() #
        self._report.Rdependencies = []  # just to speed up report.

        self.results = BooleanResults()  # for time T1
        self.results2 = BooleanResults()  # for time T2

        if config is None:
            self.config.General.pknmodel.value = self.pknmodel.filename
            self.config.General.data.value = self.data.filename
            p = ParamsGA2()
            #p.name = 'GA2'  # same parameters as GA but need to change the name
            self.config.add_section(p)

        self._called = []

    # !! should be same default values as in paramsGA
    @params_to_update
    def optimise(self,  NAFac=1, pmutation=0.5, selpress=1.2, popsize=50,
                 reltol=0.1, elitism=5, maxtime=60, sizefactor=0.0001,
                 time_index_1=1, maxgens=500, maxstallgens=100, ga_verbose=True):
        """Perform the optimisation and save results


        * Results are stored in :attr:`results`
        * Models with the tolerance are stored in :attr:`results.models`

        Parameters are those of a Genetic Algorithm used to perform
        the analysis.

        If you run again, it uses the previous best bitstirng.
        Set self.session.best_bitstring = None to start from the
        full network.
        """
        self.logging.info("Running the optimisation. Can take a very long"
                          "time. To see the progression, set verboseR "
                          "attribute to True")
        # update config GA section with user parameters
        self._update_config('GA', self.optimise.actual_kwargs)

        # keep track of the GA parameters, which may have been update above
        gad = self.config.GA.as_dict()


        # to be retrieved inside Python code
        """best_bitstring = res$bString
        best_score = res$bScore
        all_scores = res$stringsTolScores
        all_bitstrings = res$stringsTol
        reactions = model$reacID
        results = as.data.frame(res$results)
        stimuli = as.data.frame(cnolist@stimuli)
        inhibitors = as.data.frame(cnolist@inhibitors)
        species = colnames(cnolist@signals[[1]])
        optim1 = T
        """
        results.columns = [x.strip() for x in results.columns]

        columns_int = ['Generation', 'Stall_Generation']
        columns_float = ['Best_score', 'Avg_Score_Gen', 'Best_score_Gen', 'Iter_time']
        results[columns_int] = results[columns_int].astype(int)
        results[columns_float] = results[columns_float].astype(float)

        # cnograph created automatically from the reactions
        try:
            N = len(self.session.best_bitstring)
            all_bs = self.session.all_bitstrings
            df = pd.DataFrame(all_bs, columns=self.reactions_r)
            models = BooleanModels(df)
            # flatten to handle exhaustive
            import numpy as np
            models.scores = np.array(list(pylab.flatten(self.session.all_scores)))
            try:
                models.cnograph.midas = self.data.copy()
            except Exception as err:
                # does not work with ExtLiverPCB
                # CNOError: 'The cues IFNg was found in the MIDAS file but is not present in the model. Change your model or MIDAS file.'
                print("something wrong in the copying of the midas into cnograph(models)")
                print(err.message)
        except:
            N = len(self.session.best_bitstring)
            all_bs = self.session.all_bitstrings
            if N == len(self.reactions_r) and len(self.session.all_bitstrings)>0:
                df = pd.DataFrame([self.session.all_bitstrings],
                              columns=self.reactions_r)
                models = BooleanModels(df)
                models.scores = easydev.to_list(self.session.all_scores)
                self._models = models
            else:
                df = pd.DataFrame(columns=self.reactions_r)
                models = BooleanModels(df)
                models.scores = easydev.to_list(self.session.all_scores)
                self._models = models

            models.cnograph.midas = self.data.copy()

        results = {
                'best_score': self.session.best_score,
                'best_bitstring': self.session.best_bitstring,
                'all_scores': self.session.all_scores,
                'all_bitstrings': self.session.all_bitstrings,
                'reactions': self._reac_cnor2cno(self.reactions_r),
                'sim_results': self.session.sim_results,  # contains mse and sim at t0,t1,
                'results': results,
                'models': models,
                'stimuli': self.session.stimuli.copy(),
                'inhibitors': self.session.inhibitors.copy(),
                'species': self.session.species,
        }

        results['pkn'] = self.pknmodel
        results['midas'] = self.data

        self.results.results = results
        self.results.models = models

        self._called.append('optimise')

    def optimise2(self, NAFac=1, pmutation=0.5, selpress=1.2, popsize=50,
                 reltol=0.1, elitism=5, maxtime=60, sizefactor=0.0001,
                maxgens=500, maxstallgens=100, ga_verbose=True,
                  time_index_2=3):
        """

        :param NAFac:
        :param pmutation:
        :param selpress:
        :param popsize:
        :param reltol:
        :param elitism:
        :param maxtime:
        :param sizefactor:
        :param maxgens:
        :param maxstallgens:
        :param ga_verbose:
        :param time_index_2: indices as R expects it. 3, means third time points.
            time_index 1 is time zero. time_index 2 is T1, time_index 3 is T2
            TODO decrement by 1 to be Python like syntax
        :return:
        """
        raise NotImplementedError

    def plot_errors(self, close=False, show=False):
        """Plots RMSE between the data and simulated data

        The simulated data uses the best bitstring obtained
        after calling :meth:`optimise` and :meth:`optimise2`.

        If :meth:`optimise2` is called, 2 time points are used
        to plot the errors. Otherwise only 1 is used.

        """
        # todo show parameter
        assert hasattr(self.results, "_results")

        if hasattr(self.results2, "_results"):
            return self._plot_errors2(close=close, show=show)
        else:
            return self._plot_errors1(close=close, show=show)

    def _plot_errors1(self, close=False, show=False):
        raise NotImplementedError
        #midas.plot(mode="mse")

    def _plot_errors2(self, close=False, show=False):
        raise NotImplementedError

    def simulate(self, bs=None, compression=True, expansion=True):
        """Return the score of the objective function

        :param list bs: a bitstring. Must be a list of zeros and ones
            with order identical to :meth:`reactions_r`, which 
            is populated once :meth:`optimise` is called.

        """
        raise NotImplementedError

    def simulate2(self, bs1=None, bs2=None, compression=True, expansion=True):
        """Return the score of the objective function using 2 time points

        :param list bs1: a bitstring. Must be a list of zeros and ones
            with order identical to :meth:`reactions_r`, which 
            is populated once :meth:`optimise` is called.
        :param list bs2: the second bitstring. Must be a list of zeros and ones
            with length equal to the number of zeros in the first btstring.

        """
        raise NotImplementedError

    def _get_models(self):
        return self.results.models
    models = property(_get_models)

    def create_report_images(self):
        # ust a simple example of settinh the uniprot url
        # should be part of cellnopt.core
        #for node in self._pknmodel.nodes():
        #    self._pknmodel.node[node]['URL'] = "http://www.uniprot.org/uniprot/?query=Ras&sort=score"
        if "optimise" not in self._called:
            raise ValueError("You did not call any optimisation so far. Call optimise()")
        self.logging.info("Creating figures")
        self._pknmodel.plot(filename=self._report._make_filename("pknmodel.svg"),
                            show=False)
        self._pknmodel.plot(filename=self._report._make_filename("pknmodel.png"),
                            show=False)

        model = self.cnograph.copy()
        model.preprocessing()
        model.plot(filename=self._report._make_filename("expmodel.png"), show=False)

        self.plot_optimised_model(filename=self._report._make_filename("optimised_model.png"),
                                  show=False)

        self.logging.info("Creating mapback model")
        self.plot_mapback_model()
        self._report.savefig("optimised_model_mapback.png")

        self.logging.info("Creating error figure")
        self.plot_errors(show=False)
        self._report.savefig("Errors.png")

        self.logging.info("Creating midas figure")
        self.midas.plot()
        self._report.savefig("midas.png")
        pylab.close()


        self.plot_fitness(show=False, save=True)
        if "optimise2" in self._called:
            self.plot_fitness2(show=False, save=True)

    def plot_fitness2(self, show=True, save=False):
        self.results2.plot_fit()

        if save is True:
            self._report.savefig("fitness2.png")

        if show is False:
            pylab.close()

    def create_report(self):
        """Creates a full report with figures and HTML page"""

        self._create_report_header()


        txt = """<pre class="literal-block">\n"""
        # txt += "\n".join([x for x in self._script_optim.split("\n") if "write.csv" not in x])

        txt = """<a href="http://www.cellnopt.org/">
            <object data="pknmodel.svg" type="image/svg+xml" >
            <span>Your browser doesn't support SVG images</span>

            </object></a>"""

        # <img src="pknmodel.png">

        self._report.add_section(txt, "PKN graph", [("http://pythonhosted.org//cno/references.html#module-cno.io.cnograph", "cnograph")])

        txt = """<a class="reference external image-reference" href="scripts/exercice_3.py">
<img alt="MIDAS" class="align-left" src="midas.png" /></a>"""
        self._report.add_section(txt, "Data", [("http://pythonhosted.org//cno/references.html#module-cno.io.midas.xmidas", "XMIDAS class")])

        self._report.add_section('<img class="figure" src="expmodel.png">',
            "Expanded before optimisation")
        self._report.add_section('<img class="figure" src="optimised_model.png">',
            "Optimised model")

        from cno.core.report import HTMLTable
        t =  HTMLTable(self._get_model_as_df(), 'bitstrings')
        self._report.add_section(t.to_html(index=True), "Reactions on")

        self._report.add_section('<img class="figure" src="optimised_model_mapback.png">', "Optimised mapback model")

        if "optimise2" in self._called:
            self._report.add_section('<img class="figure" src="fitness.png"><img class="figure" src="fitness2.png">',
            "Fitness (T1 and T2)")

        else:
            self._report.add_section('<img class="figure" src="fitness.png">',
            "Fitness")

        self._report.add_section('<img class="figure" src="Errors.png">',
            "Errors")

        # self._report.add_section(self.get_html_reproduce(), "Reproducibility")
        fh = open(self._report.report_directory + os.sep + "rerun.py", 'w')
        fh.write("from cellnopt import CNObool, cnodata*\n")
        fh.write("CNObool(config=config.ini)\n")
        fh.write("c.optimise()\n")
        fh.write("c.onweb()\n")
        fh.close()

        # some stats
        stats = self._get_stats()
        txt = "<table>\n"
        for k in sorted(stats.keys()):
            txt += "<tr><td>%s:</td><td>%s</td></tr>\n" % (k, stats[k])
        txt += "</table>\n"
        #txt += """<img id="img" onclick='changeImage();' src="fit_over_time.png">\n"""
        self._report.add_section(txt, "stats")

        # dependencies
        self._report.add_section('<img src="dependencies.svg">', "Dependencies")
        self._report.write("index.html")

    def _get_model_as_df(self):
        assert "optimise" in self._called

        bs1 = self.results.results.best_bitstring
        N = len(bs1)

        reactions1 = list(self.results.results.reactions)
        df = pd.DataFrame({'t1': bs1, 't2': [None]*N,
                           'reactions(cno)': reactions1,
                           'reactions(CellNOptR)': self.reactions_r},
                          index=reactions1)

        if "optimise2" in self._called:
            for x, y in zip(self.results2.results.best_bitstring,
                            self.results2.results.reactions):
                df = df.set_value(y, 't2', x)
        return df

    def _get_stats(self):
        res = {}
        # res['Computation time'] = self.total_time
        try:
            res['Best Score'] = self.results.results.best_score
            res['Total time'] = self.results.results.sults.Iter_time.sum()
            res['Best Score T2'] = self.results2.results.best_score
            res['Total time T2'] = self.results2.results.results.Iter_time.sum()
        except:
            pass
        return res


def standalone(args=None):
    """This function is used by the standalone application called cellnopt_boolean

    ::

        cellnopt_boolean --help

    """
    if args is None:
        args = sys.argv[:]

    from cno.core.standalone import Standalone
    user_options = OptionsBoolean()
    stander = Standalone(args, user_options)

    # just an alias
    options = stander.options

    if options.onweb is True or options.report is True:
        trainer = CNObool(options.pknmodel, options.data, verbose=options.verbose,
            verboseR=options.verboseR, config=options.config_file, use_cnodata=options.cnodata)
        trainer.preprocessing() # should be called
    else:
        stander.help()

    trainer.optimise(**stander.user_options.config.GA.as_dict())

    if stander.user_options.config.GA2.time_index_2.value != -1:  ## need to use the time_index_2
        trainer.optimise2(**stander.user_options.config.GA2.as_dict())

    # required to call report() afterwards
    stander.trainer = trainer
    stander.report()


class OptionsBoolean(OptionsBase):
    def __init__(self):
        prog = "cno_boolean_steady"
        version = prog + " v1.0 (Thomas Cokelaer @2014)"
        super(OptionsBoolean, self).__init__(version=version, prog=prog)
        self.add_section(ParamsGA())
        self.add_section(ParamsGA2())


if __name__ == "__main__":
    """Used by setup.py as an entry point to :func:`standalone`"""
    standalone(sys.argv)
