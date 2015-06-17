import numpy as np


class Optimisation(object):
    """


    """
    def __init__(self):

        self.parameters = []
        self.scores = []

        self._best_score = None
        self._best_parameters = None

        self.guess = None

    def _get_best_score(self):
        return self._best_score
    best_score = property(_get_best_score)

    def _get_best_parameters(self):
        return self._best_parameters
    best_parameters = property(_get_best_parameters)

    def evaluator(self, candidate, args):
        """Return fitness"""
        pass

    def generator(self):
        """Return a candidate"""
        pass

    def observer(self):
        pass