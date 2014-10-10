import os
from . import __path__

name = 'ToyPB_True'

__all__ = ['description', 'model_filename', 'data_filename', 'name']

model_filename = __path__[0] + os.sep + "PKN-{0}.sif".format(name)
data_filename = __path__[0] + os.sep + "MD-{0}.csv".format(name)
description = open(__path__[0] + os.sep + "README.rst").read()

__doc__ += description


def plot():
    from cellnopt.core import CNOGraph
    CNOGraph(model_filename, data_filename).plot()
