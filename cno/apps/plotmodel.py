# -*- python -*-
#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/cno
#
##############################################################################
import os
from optparse import  OptionParser, OptionGroup
import argparse

from cno import CNOGraph, cnodata

__all__ = ["plotmodel"]


def plotmodel(args=None):
    """This function is used by the standalone application called cellnopt_plotmodel

    ::

        cellnopt_plotmodel --help

    """
    import sys
    if args == None:
        args=sys.argv[:]
    user_options = OptionCNO(prog="cellnopt_plotmodel")
    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    else:
        options = user_options.parse_args(args[1:])

    if options.cnodata is True:
        cno = CNOGraph(cnodata(options.model[0]), cnodata(options.data[0]), 
                verbose=options.verbose)
    else:
        cno = CNOGraph(options.model[0], options.data[0], verbose=options.verbose)

    compression = True
    expansion = True
    cutNONC = True

    if options.no_expansion == True: 
        expansion = False

    if options.no_compression == True:
        compression = False

    if options.no_cutnonc == True:
        cutNONC = False

    cno.preprocessing(compression=compression, expansion=expansion, cutnonc=cutNONC)

    try:
        cno.plot(viewer=options.viewer)
    except:
        print("use --show to show the image and --viewer to specify the viewer application")


class OptionCNO(argparse.ArgumentParser):

    def  __init__(self, version="1.0", prog=None):
        usage = """usage: python %s --data ToyModelMMB.csv --model ToyModelMMB.sif""" % prog
        super(OptionCNO, self).__init__(usage=usage, version=version, prog=prog)
        self.add_input_options()

    def add_input_options(self):
        """The input oiptions.

        Default is None. Keep it that way because otherwise, the contents of
        the ini file is overwritten in :class:`apps.Apps`.
        """

        group = self.add_argument_group("Inputs", 
                    """This section allows to provide path and file names of the input data.
                    If path is provided, it will be used to prefix the midas and sif filenames.
                        --path /usr/share/data --sif test.sif --midas test.csv
                    means that the sif file is located in /usr/share/data.
                    """)

        group.add_argument("--model", dest='model', 
                         default=None, type=str, nargs=1, # at most 1 model expected
                         help="Path to model (SIF format).")
        group.add_argument("--data", dest='data', nargs="+", # a least one data file expected
                         default=None, type=str,
                         help="Name of the data files")
        group.add_argument("--verbose", dest='verbose', 
                         action="store_true", 
                         help="verbose option.")
        #group.add_argument("--overwrite", dest='overwrite', 
        #                 action="store_true", 
        #                 help="overwrit existing files.")
        group.add_argument("--no-expansion", dest='no_expansion', 
                         action="store_true", 
                         help="do not expand and gates.")
        group.add_argument("--no-compression", dest='no_compression', 
                         action="store_true", 
                         help="verbose option.")
        group.add_argument("--no-cutnonc", dest='no_cutnonc', 
                         action="store_true", 
                         help="verbose option.")
        group.add_argument("--viewer", dest='viewer', default="browse",
                         action="store_true", 
                         help="verbose option.")
        group.add_argument("--use-cnodata", dest='cnodata',
                         action="store_true", 
                         help="fetch data from cno repository.")
        group.add_argument("--save-dot", dest='save_dot', 
                         action="store_false",
                         help="Save format in dot format.")



if __name__ == "__main__":
    """Used by setup.py as an entry point to :func:`standalone`

    Do not remove !!
    """
    import sys
    plotmodel(sys.argv)



