"""Main entry point to cellnopt


::

    import cellnopt as cno
    from cno import CNOGraph, XMIDAS, cnodata

    xm = XMIDAS(cnodata("MD-ToyPB.csv"))
    c = CNOGraph("PKN-ToyPB.sif", xm)
    c.plot()


"""
from __future__ import absolute_import
__version__ = '0.0.1'


# clashes between cellnopt and availabl cellnopt packages.
# Very unfortunate...
#import cellnopt.core as cnocore

from .feeder import Feeder
from .datasets import cnodata
from .misc import CNOError
