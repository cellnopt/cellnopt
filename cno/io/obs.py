"""Module to convert MIDAS file format into OBS format used in ASP

OBS files are used to check sign consistency within ASP. The output format is::

    specy = +
    specy = -

specy can be repeated.


"""
from cno import CNOError

import numpy as np

__all__ = ["OBS"]


class OBS(object):
    """Reads MIDAS file and save into observations


     ASP OBS format is of the form::

            A = +
            A = -


    .. warning:: not for production.Need to check meaning of + and - sign. It
        is positivness or slope from t1 to t2 ?

    """
    def __init__(self, data=None):
        """.. rubric:: constructor


        :param str filemame:


        """
        if data is not None:
            # local import on purpose to avoid cycling effects
            from cno import XMIDAS
            self._midas = XMIDAS(data)
        else:
            self._midas = None
        self.species = []
        self.signs = []

    def _get_midas(self):
        return self._midas
    def _set_midas(self, midas):
        self._midas = midas
    midas = property(_get_midas, _set_midas, doc='getter/setter')

    def set_time(self, time=None):
        """Simple conversion to OBS format

        For now, the sign is based on the value of the data.

        If data is positive, sign is +. Otherwise, it is set to -.

        THis may not be what is requested. Instead, we may want to focus
        on the fact that data increases or increases at a given time as compared
        to time zero.
        """
        if self.midas is None:
            raise CNOError("set midas attribute first")
        if time is None:
            time = self.midas.times[1] # time zero is useless
        if time not in self.midas.times:
            raise CNOError("time %s not found in the midas file. " % time + \
                             "Valid values are %s" % self.midas.times)

        self.species = []
        self.signs = []
        # ignore the data with time set to 0
        print("Creating obs file for time %s" % time)
        df = self.midas.df.query("time==@time")
        print(df)
        for col in df.columns:
            for v in df[col].dropna().values:

                v = float(v)
                self.species.append(col)
                if v >= 0:
                    self.signs.append("+")
                else:
                    self.signs.append("-")

    def save(self, filename):
        """Write nodes and signs into a OBS format

        ASP OBS format is of the form::

            A = +
            A = -

        """
        if len(self.species)>0:
            h = open(filename, "w")
            for n1,s in zip(self.species, self.signs):
                h.write("%s = %s\n" % (n1, s))
            h.close()
        else:
            print("nothing to save. Is the MIDAS file correct ?")





