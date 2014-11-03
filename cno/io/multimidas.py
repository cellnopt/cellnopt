# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 16:40:58 2014

@author: cokelaer
"""
from cno.io.midas import XMIDAS


class MultiMIDAS(object):
    """Data structure to store multiple instances of MIDAS files


    You can read a MIDAS file that contains several cell lines:
    and acces to the midas files usig their cell line name

    .. doctest::

        >>> mm = MultiMIDAS(cnodata("EGFR-ErbB_PCB2009.csv"))
        >>> mm.cellLines
        ['HepG2', 'PriHu']
        >>> mm["HepG2"].namesCues
        ['TGFa', 'MEK12', 'p38', 'PI3K', 'mTORrap', 'GSK3', 'JNK']

    where the list of cell line names is available in the :attr:`cellLines`
    attribute.

    Or you can start from an empty list and add instance later on using :meth:`addMIDAS`
    method.

    """
    def __init__(self, filename=None):
        """.. rubric:: constructor

        :param str filename: a valid MIDAS file (optional)

        """
        self._midasList = []
        self._names = []
        if filename:
            self.readMIDAS(filename)

    def addMIDAS(self, midas):
        """Add an existing MIDAS instance to the list of MIDAS instances

        .. doctest::

            >>> from cellnopt.core import *
            >>> m = MIDASReader(cnodata("MD-ToyPB.csv"))
            >>> mm = MultiMIDAS()
            >>> mm.addMIDAS(m)

        """
        if midas.celltypeName not in self._names:
            self._midasList.append(midas)
            self._names.append(midas.celltypeName)
        else:
            raise ValueError("midsa with same celltype already in the list")

    def readMIDAS(self, filename):
        """read MIDAS file and extract individual cellType/cellLine

        This function reads the MIDAS and identifies the cellLines. Then, it
        creates a MIDAS instance for each cellLines and add the MIDAS instance to the
        :attr:`_midasList`. The MIDAS file can then be retrieved using their
        cellLine name, which list is stored in :attr:`cellLines`.

        :param str filename: a valid MIDAS file containing any number of cellLines.


        """
        raise NotImplementedError
        m = XMIDAS(filename)
        self.addMIDAS(m)


    def _get_cellLines(self):
        names = [x.celltypeName for x in self._midasList]
        return names
    cellLines = property(_get_cellLines,
        doc="return names of all cell lines, which are the MIDAS instance identifier ")

    def __getitem__(self, name):
        index = self.cellLines.index(name)
        return self._midasList[index]

    def plot(self):
        """Call plot() method for each MIDAS instances in different figures

        More sophisticated plots to easily compare cellLines could be
        implemented.

        """
        for i,m in enumerate(self._midasList):
            from pylab import figure, clf
            figure(i+1)
            clf()
            m.plot()
