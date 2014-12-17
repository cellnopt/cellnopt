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
from biokit.rtools import RSession
from cno.core.params import CNOConfig
from easydev import Logging

__all__ = ["CNOBase", "CNORBase"]


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


class CNOBase(Logging ):
    """Alias to CNOGraph and common class to all simulators"""

    def __init__(self, pknmodel, data, tag=None, verbose=False, config=None):
        super(CNOBase, self).__init__(level=verbose)

        # TODO: check that files do exist and raise an error otherwise
        self._pknmodel = None
        self._data = None
        self.name = self.__class__.__name__
        if tag is not None:
            self.tag = tag
        else:
            self.tag = ""

        # DON'T MOVE those imports to prevent import cycling
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

        self.config = CNOConfig()

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

    def save_config_file(self, filename=None):
        if filename is None:
            filename = self.report_directory + os.sep + "config.ini"
        self.config.save(filename)

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


