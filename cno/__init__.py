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



# clashes between cellnopt and available cellnopt packages.
# Very unfortunate...
# import cellnopt.core as cnocore

from .feeder import Feeder

from .datasets import cnodata

from .io.reactions import Reactions, Reaction
from .io.sif import SIF
from .io.cna import CNA
from .io.midas import XMIDAS

from .misc import CNOError

from .testing import getdata

