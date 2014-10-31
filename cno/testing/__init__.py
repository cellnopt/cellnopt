#from __future__ import absolute_import

import os
from . import __path__
from cno.misc import CNOError

__all__ = ['getdata']


def getdata(filename):
    filename = __path__[0] + os.sep + filename
    if os.path.exists(filename):
        return filename
    else:
        import glob
        valid_filenames = glob.glob(__path__[0] + os.sep + '*')
        txt = "Unknown file {0}\n Valid files are:\n".format(filename)
        for this in valid_filenames:
            txt += " - {0}\n".format(this)
        raise CNOError(txt)
