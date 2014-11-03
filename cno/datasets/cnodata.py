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


def cnodata(filename=None):
    if filename is None:
        print("Valid names are:")
        print("\n".join(sorted(registers)))
        return

    msg = "Unknown filename. "
    msg += "Type cnodata() without argument to get the list"

    if filename in registers:
        tag = filename.split("-",1)[1].split('.')[0]

        mod = importlib.import_module('cno.datasets.{0}'.format(tag))
        if filename.startswith("PKN-"):
            fullpath = mod.model_filename
        elif filename.startswith("MD-"):
            fullpath = mod.data_filename
        else:
            print(msg)
        return fullpath
    else:
        print(msg)
