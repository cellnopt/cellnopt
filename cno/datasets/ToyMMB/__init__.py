import os
from cellnopt.core import XMIDAS, CNOGraph
from . import __path__


tag = 'ToyMMB'
__all__ = ['model', 'data', 'description', 'model_filename', 'data_filename']

try:
    model_filename = __path__[0] + os.sep + "PKN-{0}.sif".format(tag)
    model = CNOGraph(model_filename)

    data_filename = __path__[0] + os.sep + "MD-{0}.csv".format(tag)
    data = XMIDAS(data_filename)

    description = open(__path__[0] + os.sep + "README.rst").read()
except:
    print("Could not load {0} data set".format(tag))



