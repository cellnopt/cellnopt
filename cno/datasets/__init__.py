"""In order to add a data set automatically within cno,
add a directory.


"""
import os
import glob
from . import __path__


__all__ = ['cnodata']
# Finds the directories automatically.
# Nothing to change here below

# we will store the relevant directories here
_registered = []

# here is the list
pathname = __path__[0]
directories = [x for x in os.listdir(pathname) if os.path.isdir(os.sep.join([pathname,x]))]
assert len(directories)

# but some may not be valid
for register in directories:
    if len(glob.glob(os.sep.join([pathname, register, 'PKN-*'])))==0:
        pass
        #print("register {0} not valid (no PKN found)".format(register))
    if len(glob.glob(os.sep.join([pathname, register, 'MD-*'])))==0:
        pass
        #print("register {0} not valid (no MIDAS found)".format(register))
    else:
        _registered.append(register)


def _build_registers():
    registers = []
    for k in _registered:
        registers.append("PKN-{0}.sif".format(k))
    for k in _registered:
        filename = "PKN-{0}.xml".format(k)
        if os.path.exists(os.sep.join([pathname, k,filename])):
            registers.append("PKN-{0}.xml".format(k))
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
