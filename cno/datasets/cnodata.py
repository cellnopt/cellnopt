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

import importlib
from cno.datasets import registers

def cnodata(filename):

    if filename in registers:
        tag = filename.split("-")[1].split('.')[0]

        mod = importlib.import_module('cno.datasets.{0}'.format(tag))
        if filename.startswith("PKN"):
                fullpath = mod.model_filename
        else:
                fullpath = mod.data_filename
        return fullpath
    else:
        print("Unknown filename. Registered filenames are {0}".format(registers))
