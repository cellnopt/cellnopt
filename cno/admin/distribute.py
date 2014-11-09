# -*- python -*-
#
#  This file is part of the cellnopt software
#
#  Copyright (c) 2011-2012 - EMBL-EBI
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  https://pypi.python.org/pypi/cellnopt.admin
#  http://pythonhosted.org/cellnopt.admin
#
##############################################################################
# $Id: distribute.py 4207 2013-12-04 13:09:16Z cokelaer $
"""This module provides tools to check out a revision of a SVN directory
into a temporary directory and to create a tar ball out of that
SVN copy and copy the final tar ball in the directory where this script is
called. Sub directories such as .git are removed. The temporary directory is
also removed.

This is coded in the :class:`DistributeRPackage`  class. This module is also executable
provided with cellnopt.admin package called cellnopt_distribute.


"""
import os
import subprocess
import tempfile
from os.path import join as pj
import glob
import time
import sys

from easydev.tools import shellcmd
from easydev import Logging
include = ["*"]


__all__ = ["DistributeRPackage"]

class DistributeRPackage(object):
    """Class to ease distribution of CellNOptR packages from SVN


    Can be used for any SVN containing valid R packages by setting the
    repository URL

        >>> d = DistributeRPackage()
        >>> d.distribute()

    You can also use the executable provided in cellnopt.admin package itself::

        cellnopt_distribute --package CNORdt --revision HEAD

    equivalent to (if you have the sources)::

        python distribute.py --package CNORdt --revision 666
    
    Version of cellnopt.admin were based on private SVN but we moved to github and
    therefore this class is now related to the github repository only.

    In practice, this is more complicated than SVN:
    
    - cannot get a nice revision number and be able to compare
      between revisions
    - Cannot checkout a sub directory (e.g., CNORdt)

    So, we will therefore build up all packages in one go and unfortunately
    add the commit hash long number as a tag....although the --short option
    seems to be a solution.
    

    .. todo:: MEIGOR
    """
    #: mapping between package name and actual directory name 
    _valid_packages = {"CellNOptR":"CellNOptR",
                       "CNORdt":pj("CNOR_dt","CNORdt"),
                       "CNORode":pj("CNOR_ode","CNORode"),
                       "CNORfuzzy":"CNOR_fuzzy",
                       "CNORfeeder":"CNOR_add_links",
                       "MEIGOR": pj('essR','MEIGOR')}

    def __init__(self,  package, build_options="--no-build-vignettes"):
        """

        :param str package: name of a valid package (e.g., CellNOptR)
        :param revision: SVN revision (default is HEAD )
        :param build_options: additional build options for R (default is
            --no-build-vignettes)

        You can also change the logging level (e.g., self.logging.debugLevel="WARNING")


        """
        self.url = "https://github.com/cellnopt/CellNOptR"

        self.exclude = [".git"]
        self.package = package
        self.package_path = 'packages' + os.sep + package
        self.dtemp = None
        self.cwd = os.getcwd()
        self.build_options = build_options
        self.logging = Logging("INFO")

    def _get_version(self):
        data = open(self.dtemp + os.sep + self.package_path + os.sep + "DESCRIPTION", "r").read()
        res = [x.split(':')[1].strip() for x in data.split("\n") if x.startswith('Version')]
        return res[0]

    def _create_temp_directory(self):
        self.dtemp = tempfile.mkdtemp()

    def _checkout_git(self):
        self.logging.info("1. Gettting the source from GIT --------------")
        if self.dtemp == None:
            self._create_temp_directory()

        cmd = """
         git init %(directory)s;
         cd %(directory)s;
         git remote add -f origin %(repo)s true;
         echo "packages/%(package_name)s" >> .git/info/sparse-checkout;
         git pull origin master
        """
        print(self.package, self.url, self.dtemp)
        cmd = cmd % {'package_name':self.package, 'repo': self.url , 
                'directory': self.dtemp}
        
        self.logging.info(cmd,)
        try:
            ret = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            ret.wait()
            if ret.poll()!=0:
                raise Exception
            self.logging.info('...done')
        except Exception:
            raise Exception

    @staticmethod
    def help():
        """Return usage help message"""
        print("\nPURPOSE:"+__doc__)
        print("USAGE: python distribute.py --package valid_package_name")
        print("")
        print("Possible package names are %s ." % DistributeRPackage._valid_packages.keys())
        #sys.exit(1)

    def distribute(self):
        """This is the main method to create package distribution.

        It follows these steps:

         #. creates temp directoy
         #. svn checkout clean revision
         #. calls R CMD build

        """
        if self.dtemp == None:
            self._create_temp_directory()
        try:
            self._checkout_git()
        except Exception, e:
            self._stop()
            raise Exception(e)

        self._build_R_package()
        self._stop()

    def _stop(self):
        import shutil
        if self.dtemp:
            shutil.rmtree(self.dtemp)

    def clean(self):
        import os
        os.remove(self.package_name_rev)


    def _get_revision(self):
        self.logging.info("Getting revision")
        try:
            cmd = """git rev-parse --short HEAD"""
            self.logging.info(cmd)
            ret = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            ret.wait()
            revision = ret.stdout.read().strip()
        except Exception, e:
            revision = self.revision_user
            raise Exception(e)
        self.logging.info("This is revision %s. Making a tar ball." % revision)
        return revision
    # no need to make it public

    def _build_R_package(self):
        # first, build the R package
        self.logging.info("2. Creating the package distribution ----------------------"),
        t0 = time.time()
        cmdR = "R CMD build %s %s" % (self.dtemp+os.sep+ self.package_path, 
                self.build_options)
        shellcmd(cmdR, verbose=True)

        import glob
        package_name = self.package + "_" + self._get_version() + ".tar.gz"
        #rename the package name
        package_name_rev = self.package + "_"+self._get_version() + "_" + self._get_revision() + ".tar.gz"
        self.logging.info("3. "),
        #shellcmd("ls %s" % (self.dtemp), verbose=True)
        #shellcmd("ls", verbose=True)
        shellcmd("mv %s %s" % (package_name, package_name_rev), verbose=True)

        t1 = time.time()
        self.logging.info("Distribution built in " + str(t1-t0) + " seconds.")
        self.logging.info("File %s available in the current directory " % package_name_rev)
        self.package_name_rev = package_name_rev



def main():
    """Main executable related to distribute

    type::

        python distribute.py --help

    or::

        cnor_distribute --help

    """
    import sys
    print("RUNNING distribute.py")
    print("===========================================")
    print("Author: T. Cokelaer")
    if len(sys.argv) != 3:
        DistributeRPackage.help()
    else:
        import tempfile
        assert sys.argv[1] == "--package", DistributeRPackage.help()
        package = sys.argv[2]
        d = DistributeRPackage(package)
        d.distribute()

if __name__ == "__main__":
    main()
