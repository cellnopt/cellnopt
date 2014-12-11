# -*- python -*-
#
#  This file is part of the CNO package
#
#  Copyright (c) 2014 - EMBL-EBI
#
#  File author(s): Thomas Cokelaer (cokelaer@ebi.ac.uk)
#
#  Distributed under the GLPv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: www.cellnopt.org
#
##############################################################################
import os
import argparse
from biokit.rtools import RSession
import easydev


__all__ = ["CNOBase", "CNORBase", "OptionBase"]


class CNORBase(object):
    def __init__(self, verboseR):
        assert verboseR in [True, False]
        self._verboseR = verboseR
        self.session = RSession(verbose=self.verboseR)

    def _get_verboseR(self):
        return self._verboseR
    def _set_verboseR(self, value):
        assert value in [True, False]
        self._verboseR = value
        self.session.dump_stdout = value
    verboseR = property(_get_verboseR, _set_verboseR)


class CNOBase(object):
    """Alias to CNOGraph and common class to all simulators"""

    def __init__(self, pknmodel, data, tag=None, verbose=False, config=None):
        # TODO: check that files do exist and raise an error otherwise
        self._pknmodel = None
        self._data = None
        self._verbose = verbose
        self.name = self.__class__.__name__
        if tag is not None:
            self.tag = tag
        else:
            self.tag = ""

        # DONT MOVE those imports to prevent import cycling
        from cno.io import CNOGraph
        from cno.io import XMIDAS
        #
        self._data = XMIDAS(data)

        self._pknmodel = CNOGraph(pknmodel)
        self._pknmodel.midas = self._data

        self._model = CNOGraph(pknmodel)
        self._model.midas = self._data
        self._model.preprocessing() #FIXME what if one decides to preprocess differently

        self._cnograph = CNOGraph(pknmodel, data)

        if config:
            self.load_config(config)
        else:
            # if use_cnodata, let us search files using cnodata if not found locally
            """if use_cnodata:
                if os.path.exists(pknmodel) == False:
                    pknmodel = cnodata(pknmodel)
                    if os.path.exists(data) == False:
                        data = cnodata(data)
                                                                                                                    if pknmodel == None or data == None:
                txt = "Input %s not a valid file" % pknmodel
                txt += "Or %s not a valid file" % data
                raise ValueError(txt)
                                                                                                                    self._pknmodel_filename = pknmodel
            self._data_filename = data"""
            # keep track of all configuration parameters
            self.init_config()

    def _get_verbose(self):
        return self._verbose
    def _set_verbose(self, verbose):
        # TODO check value is a bool
        self._verbose = verbose
    verbose = property(_get_verbose, _set_verbose)

    def _get_pknmodel(self):
        return self._pknmodel
    pknmodel = property(_get_pknmodel, doc="get the prior knowledge network")

    def _get_cnograph(self):
        return self._cnograph
    cnograph = property(_get_cnograph)

    def _get_data(self):
        return self._data
    data = property(_get_data, doc="get the data (MIDAS)")
    midas = property(_get_data, doc="get the data (MIDAS)")

    def preprocessing(self, expansion=True, compression=True, cutnonc=True,
            maxInputsPerGate=2):
        """Apply preprocessing on the PKN model"""
        # FIXME if used, the pknmodel is changed
        self._pknmodel.preprocessing(expansion, compression, cutnonc,
                maxInputsPerGate=maxInputsPerGate)

    def plot_pknmodel(self):
        """Plot the original PKN model

        .. plot::
            :include-source:
            :width: 80%

            from cellnopt.pipeline.cnobase import CNObase
            from cellnopt.data import cnodata
            o = CNOBase(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"),
                formalism="base")
            o.plot_midas(xkcd=True)   # use xkcd for fun


        .. seealso:: full documentation about MIDAS in cellnopt.core.cnograph
        """
        self._pknmodel.plot()

    def plot_model(self):
        self._model.plot()

    def plot_optimised_model(self, filename=None, show=True):
        bs = self.results.results.best_bitstring
        reactions = self.results.results.reactions
        opt = {}
        for b,r in zip(bs, reactions):
            if b == 0:
                opt[r]= 0
            else:
                opt[r] = 1
        m = self._model.copy()
        m.set_edge_attribute('opt', opt)
        m.plot(edge_attribute='opt', cmap='gray_r', show=show, filename=filename)

    def _reac_cnor2cno(self, reactions):
        # CNOR ands are encoded with +
        from cno import Reaction
        reactions = [r.replace('+', '^') for r in reactions]
        newreacs = []
        for i,reac in enumerate(reactions):
            r = Reaction(reac)
            r.sort()
            newreacs.append(r.name)
        return newreacs

    def plot_mapback_model(self):
        from cno.io import mapback
        m = mapback.MapBack(self.pknmodel, self._model)
        reactions = self.results.results.reactions
        bs = self.results.results.best_bitstring
        links2map = [r for b,r in zip(bs,reactions) if b==1]
        newreacs  = m.mapback(links2map)
        model = m.plot_mapback(newreacs)
        return model

    def plot_midas(self, xkcd = False):
        """Plot the MIDAS data

        .. seealso:: full documentation about MIDAS in :meth:`cellnopt.core.midas.plot`
        """
        if xkcd:
            from pylab import xkcd, gcf
            with xkcd():
                self.midas.plot()
                f = gcf()
                f.set_facecolor("white")
        else:
            self.midas.plot()

    def get_html_reproduce(self):
        text = """
        <p>In a shell, go to the report directory (where is contained this report)
        and either execute the file called rerun.py or typoe these commands in
        a python shell
        </p>

        <pre>
        from cellnopt.pipeline import *
        c = CNObool(config=config.ini)
        c.gaBinaryT1()
        c.report()
        </pre>

        <p>You will need the configuration file that was used in the report
        (<a href="config.ini">config.ini</a>) </p>

        <p>New results will be put in a sub directory so that your current
        report is not overwritten</p>
        """
        return text

    def load_config(self, filename):
        import easydev.config_tools
        self.config = easydev.config_tools.DynamicConfigParser()
        self.config.load_ini(filename)

    def init_config(self, force=False):
        if "config" in self.__dict__.keys():
            if force == False:
                raise ValueError("""Your config already exists and will be erased.
If you really want to re-initialise the config attribute, use force=True
parameter.""" )

        self.config = easydev.config_tools.DynamicConfigParser()
        self.config.add_section("General")
        self.config.add_option("General", "pknmodel", self._pknmodel.filename)
        self.config.add_option("General", "data", self._data.filename)
        self.config.add_option("General", "formalism", self.__class__.__name__)
        #self.config.add_option("General", "use_cnodata", self._use_cnodata)
        self.config.add_option("General", "tag", self.tag)
        #self.config.add_option("General", "overwrite_report", self._overwrite_report)
        #self.config.add_option("General", "Rexecutable", self.Rexecutable)
        self.config.add_option("General", "verbose", self.verbose)
        #self.config.add_option("General", "report_directory", self.report_directory)

    def save_config_file(self, filename=None):
        if filename==None:
            filename = self.report_directory + os.sep + "config.ini"
        self.config.save(filename)


class OptionBase(argparse.ArgumentParser):

    def  __init__(self, version="1.0", prog=None, usage=""):
        super(OptionBase, self).__init__(usage=usage, version=version, prog=prog)
        self.add_general_options()

    def add_general_options(self):
        """The input oiptions.

        Default is None. Keep it that way because otherwise, the contents of
        the ini file is overwritten in :class:`apps.Apps`.
        """

        group = self.add_argument_group("General",
                    """This section allows to provide path and file names of the input data.
                    If path is provided, it will be used to prefix the midas and sif filenames.
                        --path /usr/share/data --sif test.sif --midas test.csv
                    means that the sif file is located in /usr/share/data.
                    """)

        group.add_argument("--model", dest='model',
                         default=None, type=str, # at most 1 model expected
                         help="Path to model (SIF format).")
        group.add_argument("--data", dest='data', # a least one data file expected
                         default=None, type=str,
                         help="Name of the data files")
        group.add_argument("--config", dest='config', # a least one data file expected
                         default=None, type=str,
                         help="Name of configuration file")
        group.add_argument("--verbose", dest='verbose',
                         action="store_true",
                         help="verbose option.")
        group.add_argument("--no-expansion", dest='no_expansion',
                         action="store_true",
                         help="do not expand and gates.")
        group.add_argument("--no-compression", dest='no_compression',
                         action="store_true",
                         help="compression option.")
        group.add_argument("--no-cutnonc", dest='no_cutnonc',
                         action="store_true",
                         help="nonc option.")
        group.add_argument("--report", dest='report',
                         action="store_true",
                         help="report option.")
        group.add_argument("--use-cnodata", dest='use_cnodata',
                         action="store_true",
                         help="user cnodata to fetch file from cellnopt.data (local or web version).")
        group.add_argument("--overwrite-report", dest='overwrite_report',
                         action="store_false",
                         help="overwrite existing directory.")
        group.add_argument("--tag", dest='tag', # a least one data file expected
                         default=None, type=str,
                         help="tag to append to the report directory name")
        group.add_argument("--config-file", dest='config', # a least one data file expected
                         default=None, type=str,
                         help="a configuration file")


