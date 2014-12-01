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
from .boolean import *

# FEEDER package
from .feeder import Feeder

# SIMIULATOR
try:
    from .boolean import CNORbool

except Exception as err:
    print(err.message)
    print("Issue in boolean package. Please report the issue to github/cellnopt/cellnopt")


#MINLP
from .milp import *


#ODE
from .ode import CNORode

#fuzzy
from .fuzzy import CNORfuzzy

