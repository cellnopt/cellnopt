

_registered = ['ToyMMB', 'ToyPCB', 'ToyPB', 'ToyPB_True', 'ExtLiverPCB',
'LiverDREAM']

def _build_registers():
    registers = []
    for k in _registered:
        registers.append("PKN-{0}.sif".format(k))
    for k in _registered:
        registers.append("MD-{0}.csv".format(k))
    registers = sorted(registers)
    return registers



registers = _build_registers()
from .cnodata import cnodata
