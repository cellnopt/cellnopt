
__all__ = ["CNOBase"] 

class CNOBase(object):
    """Alias to CNOGraph and common class to all simulations"""
    def __init__(self, pknmodel, data, verbose=False):

        self._pknmodel = None
        self._data = None
        self._verbose = verbose

        # prevent import cycling
        from cno.io import CNOGraph
        from cno.io import XMIDAS
        self._data = XMIDAS(data)
        self._pknmodel = CNOGraph(pknmodel)
        self._pknmodel.midas = self._data
        
    def _get_verbose(self):
        return self._verbose
    def _set_verbose(self, value):
        # TODO check value is a bool
        self.verbose = verbose
    verbose = property(_get_verbose, _set_verbose)

    def _get_model(self):
        return self._pknmodel
    pknmodel = property(_get_model, doc="get the prior knowledge network")

    def _get_data(self):
        return self._data
    data = property(_get_data, doc="get the data (MIDAS)")

    def preprocessing(self, expansion=True, compression=True, cutnonc=True, 
            maxInputsPerGate=2):
        """Apply preprocessing on the PKN model"""
        self._pknmodel.preprocessing(expansion, compression, cutnonc, 
                maxInputsPerGate=maxInputsPerGate)

