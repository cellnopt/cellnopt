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
#  website: http://github.com/cellnopt/cellnopt
#
##############################################################################
import os

from cno.core.params import CNOConfig
from cno.io.reactions import Reaction
from cno.core.params import CNOConfigParser
from cno.datasets import cnodata

from biokit.rtools import RSession
from easydev import Logging

import pylab


__all__ = ["CNOBase", "CNORBase"]


class CNORBase(object):
    """A base class to handle the R session and its verbosity"""
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


class CNOBase(Logging):
    """Abstract Base Class common class to all formalisms"""

    def __init__(self, pknmodel, data, tag=None, verbose=False, use_cnodata=False, config=None):
        super(CNOBase, self).__init__(level=verbose)

        # DON'T MOVE those imports to prevent import cycling
        from cno.io import CNOGraph
        from cno.io import XMIDAS

        # Default PKN is not preprocessed by default
        self._expansion = False
        self._compression = False
        self._cutnonc = False
        self._max_inputs_per_gate = 3

        # Now reads the config file if any.
        # Do we want to call preprocessing()
        # TODO 
        self.config = CNOConfig()
        if config is not None:
            if isinstance(config, str):
                self.config.read(config)
            else:
                raise TypeError

        self.config.General.cnodata.value = use_cnodata

        # TODO: check that files do exist and raise an error otherwise
        self._pknmodel = None
        self._data = None
        self.name = self.__class__.__name__
        if tag is not None:
            self.tag = tag
        else:
            self.tag = ""

        # data and model may not be found. If cnodata provided, search in the package
        if self.config.General.cnodata.value is True:
            data = cnodata(data)
            pknmodel = cnodata(pknmodel)
            
        self._data = XMIDAS(data)
        self._pknmodel = CNOGraph(pknmodel)
        self._pknmodel.midas = self._data
        self._model = self._pknmodel.copy()

        # why another copy ?
        self._cnograph = self._pknmodel.copy()


        self._report_directory = None
        self.report_directory = 'report'

    def _get_report_directory(self):
        return self._report_directory
    def _set_report_directory(self, value):
        self._report_directory = value
    report_directory = property(_get_report_directory, _set_report_directory)

    def _update_config(self, section, kwargs):
        for k, v in kwargs.items():
            try:
                this_section = getattr(self.config, section)
                option = getattr(this_section, k)
                option.value = v
            except:
                # additional options should be ignored.
                pass

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
        self._model = self._pknmodel.copy()
        self._model.preprocessing(expansion, compression, cutnonc,
                maxInputsPerGate=maxInputsPerGate)
        self._expansion = expansion
        self._compression = compression
        self._cutnonc = cutnonc
        self._max_inputs_per_gate = maxInputsPerGate

    def plot_pknmodel(self):
        """Plot the original PKN model


        .. seealso:: full documentation about MIDAS in cellnopt.core.cnograph
        """
        self._pknmodel.plot()

    def plot_model(self):
        self._model.plot()

    def plot_optimised_model(self, filename=None, show=True,
            show_pruned_edges=True):
        try:
            bs = self.results.results.best_bitstring
            reactions = self.results.results.reactions
        except:
            bs = self.best_bitstring
            reactions = self.reactions2parameters(bs)
        opt = {}
        for b,r in zip(bs, reactions):
            if "!" in r: #inhibitory
                if b == 0:
                    if show_pruned_edges:
                        opt[r]= 'pink'
                    else:
                        opt[r]= 'white'
                else:
                    opt[r] = 'red'
            else:
                if b == 0:
                    if show_pruned_edges:
                        opt[r]= 'lightgray'
                    else:
                        opt[r]= 'white'
                else:
                    opt[r] = 'black'


        m = self._model.copy()
        m.set_edge_attribute('color', opt)
        #m.plot(edge_attribute='opt', cmap='gray_r', show=show, filename=filename)
        m.plot(show=show, filename=filename)
        return m

    def _reac_cnor2cno(self, reactions):
        # CNOR ands are encoded with +

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

    def save_config_file(self, filename=None):
        if filename is None:
            filename = self.report_directory + os.sep + "config.ini"
        self.config.save(filename)

    def report(self, show=True):
        self.onweb(show=show)

    def onweb(self, show=True):
        self.create_report()
        try:
            self.create_report_images()
        except:
            pass
        if show:
            self._report.show()

    def create_report(self):
        raise NotImplementedError

    def create_report_images(self):
        raise NotImplementedError

    def _create_report_header(self):
        self._report._init_report()
        self._report.directory = self._report.report_directory
        # Save filenames and report in a section
        fname = self._report.directory + os.sep + "PKN-pipeline.sif"

        self.config.save(self._report.directory + os.sep + 'config.ini')
        self.cnograph.to_sif(fname)
        fname = self._report.directory + os.sep + "MD-pipeline.csv"
        self.midas.to_midas(fname)
        txt = '<ul><li><a href="PKN-pipeline.sif">input model (PKN)</a></li>'
        txt += '<li><a href="MD-pipeline.csv">input data (MIDAS)</a></li>'
        txt += '<li><a href="config.ini">Config file</a></li>'
        txt += '<li><a href="rerun.py">Script</a></li></ul>'
        txt += "<bold>some basic stats about the pkn and data e.g. number of species ? or in the pkn section?</bold>"
        self._report.add_section(txt, "Input data files")

        self._report.add_section(
        """
         <div class="section" id="Script_used">
         <object height=120 width=300 type='text/x-scriptlet' border=1
         data="description.html"></object>
         </div>""", "Description")

    def plot_fitness(self, show=True, save=False):
        # TODO show and save parameters
        self.results.plot_fit()

        if save is True:
            self._report.savefig("fitness.png")

        if show is False:
            pylab.close()
