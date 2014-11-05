"""Main entry point to cellnopt


::

    import cellnopt as cno
    from cno import CNOGraph, XMIDAS, cnodata

    xm = XMIDAS(cnodata("MD-ToyPB.csv"))
    c = CNOGraph("PKN-ToyPB.sif", xm)
    c.plot()


"""
from __future__ import absolute_import

try:
    import pkg_resources
    version = pkg_resources.require("cno")[0].version
    __version__ = version
except:
    version = "undefined"



# DATASETS package
from .datasets import cnodata

# IO package

from .io.reactions import Reactions, Reaction
from .io.sif import SIF
from .io.cna import CNA
from .io.midas import XMIDAS
from .io.cnograph import CNOGraph

# MISC package
from .misc import CNOError

# DATA for TESTING package
from .testing import getdata


# CORE package
from .core import *
 
# FEEDER package
from .feeder import Feeder

# SIMIULATOR
from .boolean import CNORbool
