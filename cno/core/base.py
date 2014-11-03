from cellnopt.core import CNOGraph
from cellnopt.core import XMIDAS


class CNOBase(object):


    def __init__(self, model, data, verbose=False):

        self._model = None
        self._data = None
        self._verbose = verbose


        self._model = CNOGraph(model)
        self._data = XMIDAS(data)
        


    def _get_verbose(self):
        return self._verbose
    def _set_verbose(self, value):
        # TODO check value is a bool
        self.verbose = verbose
    verbose = property(_get_verbose, _set_verbose)

    def _get_model(self):
        return self._model
    model = property(_get_model)

    def _get_data(self):
        return self._data
    data = property(_get_data)
