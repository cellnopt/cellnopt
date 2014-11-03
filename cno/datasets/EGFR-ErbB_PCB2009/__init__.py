# -*- python -*-
#
#  This file is part of cellnopt software
#
#  Copyright (c) 2014 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: http://github.com/cellnopt
#
##############################################################################
""""""

import os
from . import __path__

name = 'EGFR-ErbB_PCB2009'

__all__ = ['description', 'model_filename', 'data_filename', 'name']

model_filename = __path__[0] + os.sep + "PKN-{0}.sif".format(name)
data_filename = __path__[0] + os.sep + "MD-{0}.csv".format(name)
description = open(__path__[0] + os.sep + "README.rst").read()

__doc__ += description


def plot():
    from cellnopt.core import CNOGraph
    CNOGraph(model_filename, data_filename).plot()

