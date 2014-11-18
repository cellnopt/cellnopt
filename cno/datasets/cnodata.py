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
    """Return the full path name of a registered model or data set

    A register directory is one contained in cno.datasets with a PKN
    and MIDAS file and __init__.py

    If you do not know the registered directories, type cnodata without
    parameter::

        cnodata()

    else, it returns the full path of an existing file::

        cnodata("PKN-ToyPB.sif")


    """
    if filename is None:
        print("Valid names are:")
        print("\n".join(sorted(registers)))
        return

    msg = "Unknown filename. "
    msg += "Type cnodata() without argument to get the list"

    if filename in registers:
        tag = filename.split("-", 1)[1].split('.')[0]

        mod = importlib.import_module('cno.datasets.{0}'.format(tag))
        if filename.startswith("PKN-") and filename.endswith('.sif'):
            fullpath = mod.model_filename
        elif filename.startswith("PKN-") and filename.endswith('.xml'):
            fullpath = mod.model_filename.replace('.sif', '.xml')
        elif filename.startswith("MD-"):
            fullpath = mod.data_filename
        else:
            print(msg)
        return fullpath
    elif "*" in filename:
        tags = [x for x in filename.split("*") if x]
        for register in registers:
            found = True
            for tag in tags:
                if tag not in register:
                    found = False
            if found:
                print(register)
    else:
        print(msg)



