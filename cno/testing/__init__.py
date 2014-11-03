#from __future__ import absolute_import

import os
import glob

from . import __path__

from cno.misc import CNOError

__all__ = ['getdata', 'getdata_metabolites']


def getdata(filename=None, pattern=None, verbose=True, raise_error=True):
    """Return the fullpath of a filename to be found in the test data set

    cellnopt contains a directory called cno/testing that provides
    data sets used within the test suite, which can also be multi-purpose within
    notebooks or demonstration.

    This function can be used without argument to figure out the
    list of available files. Otherwise, it tries to return the fullpath. If not
    found, a CNOError is returned.

    If filename is incorrect, the list of valid names is printed. To prevent
    that behaviour, you can set the **verbose** parameter to False.

    if pattern provided, raise_error is sset to False and list of
    filenames that contain the pattern is returned
    """
    if filename:
        filename = __path__[0] + os.sep + filename
    if filename and os.path.exists(filename):
        return filename
    else:
        import glob
        valid_filenames = glob.glob(__path__[0] + os.sep + '*')
        txt = "Unknown file {0}\n Valid files are:\n".format(filename)
        for this in valid_filenames:
            if this.endswith("pyc") or this.endswith('.py'):
                continue
            if os.path.isdir(this):
                continue
            this_end = this.split("cno" + os.sep + "testing")[1][1:]
            txt += " - {0}\n".format(this_end)

        metabolites = getdata_metabolites()
        for this in metabolites:
            this_end = this.split("cno" + os.sep + "testing")[1][1:]
            txt += " - {0}\n".format(this_end)

        print(txt)
        if pattern:
            raise_error = False
            valid_filenames = [filename for filename in valid_filenames
            if pattern.replace("*", "") in filename]
        if raise_error:
            raise CNOError("Incorrect filename. See list of valid filenames above.")
        else:
            return valid_filenames


def getdata_metabolites():
    metabolites = glob.glob(__path__[0] + os.sep
                + 'metabolites_samples' + os.sep + '*')
    return metabolites


