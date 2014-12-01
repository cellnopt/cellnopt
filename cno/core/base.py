# -*- python -*-
#
#  This file is part of the cellnopt package
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
__all__ = ["CNOBase"]


class CNOBase(object):
    """Alias to CNOGraph and common class to all simulators"""
    def __init__(self, pknmodel, data, verbose=False):
        # TODO: check that files do exist and raise an error otherwise
        self._pknmodel = None
        self._data = None
        self._verbose = verbose

        # DONT MOVE those imports to prevent import cycling
        from cno.io import CNOGraph
        from cno.io import XMIDAS
        #
        self._data = XMIDAS(data)
        self._pknmodel = CNOGraph(pknmodel)
        self._pknmodel.midas = self._data

        self._cnograph = CNOGraph(pknmodel, data)

    def _get_verbose(self):
        return self._verbose
    def _set_verbose(self, verbose):
        # TODO check value is a bool
        self._verbose = verbose
    verbose = property(_get_verbose, _set_verbose)

    def _get_model(self):
        return self._pknmodel
    pknmodel = property(_get_model, doc="get the prior knowledge network")

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



