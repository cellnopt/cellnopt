"""In order to add a data set automatically within cno, 
add a directory.


"""
import os
import glob
from . import __path__
# Finds the directories automatically. 
# Nothing to change here below

# we will store the relevant directories here
_registered = []

# here is the list
directories = [x for x in os.listdir(__path__[0]) if os.path.isdir(x)]

# but some may not be valid
for register in directories:
    if len(glob.glob(register + os.sep + 'PKN-*'))==0:
        print("register {0} not valid (no PKN found)".format(register))
    elif len(glob.glob(register + os.sep + 'MD-*'))==0:
        print("register {0} not valid (no MIDAS found)".format(register))
    else:
        _registered.append(register)
            

def _build_registers():
    registers = []
    for k in _registered:
        registers.append("PKN-{0}.sif".format(k))
    for k in _registered:
        registers.append("MD-{0}.csv".format(k))
    registers = sorted(registers)
    return registers

# add the PKN and MIDAS file in the list of files that can be
# retrieved.
registers = _build_registers()

# dynamic import of all directories that contain a PKN and MIDAS 
# files
for register in _registered:
    import importlib
    importlib.import_module('cno.datasets.{0}'.format(register))


from .cnodata import cnodata
