import importlib
from cno.datasets import registers

def cnodata(filename):

    if filename in registers:
        tag = filename.split("-")[1].split('.')[0]

        mod = importlib.import_module('cno.datasets.{0}'.format(tag))
        if filename.startswith("PKN"):
                fullpath = mod.model_filename
        else:
                fullpath = mod.data_filename
        return fullpath
    else:
        print("Unknown filename. Registered filenames are {0}".format(registers))
