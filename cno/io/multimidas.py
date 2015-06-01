# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 16:40:58 2014

@author: cokelaer
"""
from cno.io.midas import XMIDAS
import pylab

class MultiMIDAS(object):
    """Data structure to store multiple instances of MIDAS files


    You can read a MIDAS file that contains several cell lines:
    and acces to the midas files usig their cell line name

    .. doctest::

        >>> mm = MultiMIDAS(cnodata("MD-EGFR-ErbB_PCB2009.csv"))
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
        self._cellLines = []
        if filename is not None:
            self.readMIDAS(filename)

    def addMIDAS(self, midas):
        """Add an existing MIDAS instance to the list of MIDAS instances

        .. doctest::

            >>> from cellnopt.core import *
            >>> m = MIDASReader(cnodata("MD-ToyPB.csv"))
            >>> mm = MultiMIDAS()
            >>> mm.addMIDAS(m)

        """
        if midas.cellLines not in self._cellLines:
            self._midasList.append(midas.copy())
            self._cellLines.append(midas.cellLine)
        else:
            raise ValueError("midas with same celltype already in the list")

    #TODO rename with lower case
    #TODO allows a glob directly
    def readMIDAS(self, filename):
        """read MIDAS file and extract individual cellType/cellLine

        This function reads the MIDAS and identifies the cellLines. Then, it
        creates a MIDAS instance for each cellLines and add the MIDAS instance to the
        :attr:`_midasList`. The MIDAS file can then be retrieved using their
        cellLine name, which list is stored in :attr:`cellLines`.

        :param str filename: a valid MIDAS file or a list of filenames.

        """
        if isinstance(filename, str):
            self._read_midas(filename)
        elif isinstance(filename, list):
            for this in filename:
                self._read_midas(this)

    def _read_midas(self, filename):
        try:
            m = XMIDAS(filename)
        except:
            pass
        for cellLine in m.cellLines:
            m.cellLine = cellLine
            self.addMIDAS(m)

    def _get_cellLines(self):
        names = [x.cellLine for x in self._midasList]
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
            pylab.title('%s' % m.cellLine, fontsize=16)

    def _get_species(self):#
        if len(self.cellLines):
            cell = self.cellLines[0]
            return self[cell].df.columns[:]
        else:
            return []
    species = property(_get_species)



    def get_midas_for_one_species(self, species):
        """


        
        :param species: name of a species (i.e. column). All column that have that name
            will be merged into a new MIDAS files. new columns will be the cell line names
            and data will be from the species only.
        :return: a MIDAS file. The name of the cell line is replaced wit the species name
            and columns are the cell lines found in the MultiMIDAS structure.
        """

        import pandas as pd
        cellLines = self.cellLines[:]
        df = pd.concat([self[cell].df[species] for cell in cellLines]).unstack(0)

        from cno.io.measurements import MIDASBuilder, Measurement, Measurements


        mb = MIDASBuilder()
        # loop over the data looking at the row (experiment name and time)
        for exp, time in df.index:
            # and the column, supposibly the cell line name
            for column in df.columns:
                # fetching (1) the experiment details
                exp_details = self[cellLines[0]].experiments.ix[exp]
                # and the value
                datum = df[column].ix[exp].ix[time]
                # Create the experiment instance
                exp_instance = Measurement(column, time,
                    stimuli=exp_details.Stimuli.to_dict(),
                    inhibitors=exp_details.Inhibitors.to_dict(),
                    value=datum, cellLine=species)
                mb.add_measurements(exp_instance)
        newm = mb.xmidas
        newm.sort_experiments_by_inhibitors()
        return newm


