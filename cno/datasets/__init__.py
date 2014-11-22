"""In order to add a data set automatically within cno,
add a directory.


"""
import os
import glob
from . import __path__


__all__ = ['cnodata']
# Finds the directories automatically.
# Nothing to change here below


class Register(object):
    """Retrieve existing Network and MIDAS files in :mod:`cno.datasets`


        # finds all directories
        r = Register() # finds all directories
        # keep those that are valid in the _registered attribute
        r.filter_directories() # keep those that are valid in the _registered attribute

    """
    def __init__(self, pathname=__path__[0], debug=False):
        self.debug = debug
        self.pathname = pathname[:]

        self.directories = [x for x in os.listdir(pathname) 
                if os.path.isdir(os.sep.join([pathname,x]))]

        self.directories = sorted(self.directories)
        assert len(self.directories)

        self.filter_directories()
        self.registered = []

    def filter_directories(self):
        # but some may not be valid
        valid = []
        for register in self.directories:
            Npkn = len(glob.glob(os.sep.join([self.pathname, register, 'PKN-*'])))
            Nmidas = len(glob.glob(os.sep.join([self.pathname, register, 'MD-*'])))
            if Npkn == 0:
                if self.debug:
                    print("CNO warning (datasets): {0} directory has no PKN".format(register))
                    print("CNO warning (datasets): {0} will not be available (no PKN)".format(register))
                continue
            if Nmidas == 0:
                if self.debug:
                    print("CNO warning (datasets): {0} directory has no MIDAS".format(register))
            valid.append(register)
        self.directories = sorted(list(set(valid[:])))

    def import_all(self):
        self.registered = {}
        class _PLOT(object):
            def __init__(self, model, data):
                self.model = model
                self.data = data
            def plot(self, **kargs):
                from cno.io import CNOGraph
                CNOGraph(self.model, self.data).plot(**kargs)

        for package in self.directories:
            if self.debug:
                print("Importing %s" % package)
            import importlib
            try:
                this = importlib.import_module('cno.datasets.{0}'.format(package))
                # attach README dynamically
                this.__doc__ = open(this.__path__[0] + os.sep + 'README.rst', 'r').read()

                try:
                    metadata = this.metadata
                except:
                    metadata = {}

                if 'name' in metadata.keys():
                    name = metadata['name']
                else:
                    name = os.path.split(this.__path__[0])[1]

                if 'model' in metadata.keys():
                    # replace / by os.sep for multiplatform compat
                    model = metadata['model'].replace('/', os.sep)
                else:
                    model = "PKN-" + name + ".sif" 

                if 'data' in metadata.keys():
                    # replace / by os.sep for multiplatform compat
                    data = metadata['data'].replace('/', os.sep)
                else:
                    data = "MD-" + name + ".csv" 

                this.model = this.__path__[0] + os.sep + model
                this.data = this.__path__[0] + os.sep + data
                this.plot = _PLOT(this.model, this.data).plot

                sbmlpath = this.model.replace(".sif", ".xml")
                sbml = model.replace(".sif", ".xml")

                self.registered[model] = this.model

                # If .. is found, it means it is an existing MIDAS file
                # in another directory.                  
                if ".." not in this.data:
                    self.registered[data] = this.data
                else:
                    data = "MD-" + package + '.csv'
                    self.registered[data] = this.data

                if os.path.exists(sbmlpath):
                    self.registered[sbml]  = sbmlpath

            except Exception as err:
                print(err.message)
                print("CNO warning (datasets): could not import {0}".format(package))
                pass


# Now, we build an instance of Register that filters out directories, which do not have
# at least one PKN, and finally import dynamically all packages in the namespace
register = Register(__path__[0])
register.filter_directories()
register.import_all()
names = register.directories
# register variable is used inside cnodata

# import cnodata function in the dataset package
from .cnodata import cnodata


