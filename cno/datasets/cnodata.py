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
from . import register


def cnodata(filename=None):
    """Return the full path name of a registered model or data set

    A register directory is one contained in cno.datasets with a PKN
    and MIDAS file and __init__.py

    If you do not know the registered directories, type cnodata without
    parameter::

        cnodata()

    else, it returns the full path of an existing file::

        cnodata("PKN-ToyPB.sif")

    You can also search for a pattern::

        cnodata('*SBML*')

    """
    if filename is None:
        print("Valid names are:")
        print("\n".join(sorted(register.registered.keys())))
        return

    msg = "Unknown filename. "
    msg += "Type cnodata() without argument to get the list"

    if filename in register.registered.keys():
        return register.registered[filename]
    elif "*" in filename:
        tags = [x for x in filename.split("*") if x]
        for candidate in register.registered.keys():
            found = True
            for tag in tags:
                if tag not in candidate:
                    found = False
            if found:
                print(candidate)
    else:
        print(msg)



