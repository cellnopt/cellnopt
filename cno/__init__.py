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

# IO package. No need to expose everything.
from .io import midas, cnograph, sif # 3 most common modules.  
from .io.reactions import Reactions, Reaction
from .io.sif import SIF
from .io.eda import EDA
from .io.cna import CNA
from .io.midas import XMIDAS
from .io.cnograph import CNOGraph
from .io.xcnograph import XCNOGraph

# MISC package
from .misc import CNOError

# DATA for TESTING package
from .testing import getdata

# CORE package
from .core import *

#Boolean package
from .boolean import cnorbool, CNORbool, steady, CASPO, asp

# FEEDER package
from .feeder import Feeder

# ODE
from .ode import cnorode, CNORode

# fuzzy
from .fuzzy import cnorfuzzy, CNORfuzzy

# discrete time
from .boolean import cnordt, CNORdt

#MINLP
from .milp import *
