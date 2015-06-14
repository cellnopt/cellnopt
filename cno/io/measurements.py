# -*- python -*-
#
#  This file is part of the cinapps.tcell package
#
#  Copyright (c) 2012-2013 - EMBL-EBI
#
#  File author(s): Thomas Cokelaer (cokelaer@ebi.ac.uk)
#
#  Distributed under the GLPv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: www.cellnopt.org
#
##############################################################################
from __future__ import print_function
from __future__ import unicode_literals
import random
import copy

import numpy as np
import pandas as pd


__all__ = ["Measurement", "Measurements", "MIDASBuilder"]


class Measurement(object):
    """Data structure to store a measurement.

    Givem a list of stimuli and inhibitor, stores a measure
    at a given time.


        >>> from cno.io.midas_extra import Measurement
        >>> m = Measurement("AKT", 0, {"EGFR":1}, {"AKT":0}, 0.1)

    """
    def __init__(self, protein_name, time, stimuli, inhibitors, value,
            cellLine="undefined", units='second'):
        """

        :param str protein:
        :param float time:
        :param dict stimuli: a dictionary
        :param dict inhibitors: a dictionary
        :param float measurement: the value
        :param str cellLine: Defaults to "undefined"
        :param str units: Defaults to "second" (not yet used)

        """
        self._time = time

        self._stimuli = None
        self.stimuli = stimuli

        self._measurement = value
        self._protein_name = protein_name

        self._inhibitors = None
        self.inhibitors = inhibitors

        self._cellLine = cellLine
        self._units = units

    def _get_units(self):
        return self._units
    def _set_units(self, units):
        assert units in ["second", "hour", "minute", "day"]
        self._units = units
    units = property(_get_units, _set_units ,doc="units (second, hour, minute, day")

    def _get_cellLine(self):
        return self._cellLine
    def _set_cellLine(self, cellLine):
        assert isinstance(cellLine, str)
        self._cellLine = cellLine
    cellLine = property(_get_cellLine, _set_cellLine)

    def _get_time(self):
        return self._time
    def _set_time(self, time):
        assert isinstance(time, (int, float))
        assert time>=0
    time = property(_get_time, _set_time)

    def _get_protein_name(self):
        return self._protein_name
    def _set_protein_name(self, name):
        assert isinstance(name, str)
        self._protein_name = name
    protein_name = property(_get_protein_name, _set_protein_name)

    def _get_stimuli(self):
        return self._stimuli
    def _set_stimuli(self, stimuli):
        isinstance(stimuli, dict)
        for k,v in stimuli.items():
            if v >1 or v<0:
                raise ValueError("Value of the stimulus {} must be inside the range [0,1]".format(k))
        self._stimuli = stimuli.copy()
    stimuli = property(_get_stimuli, _set_stimuli)

    def _get_inhibitors(self):
        return self._inhibitors
    def _set_inhibitors(self, inhibitors):
        isinstance(inhibitors, dict)
        for k,v in inhibitors.items():
            if v >1 or v<0:
                raise ValueError("Value of the inhibitor {} must be inside the range [0,1]".format(k))
        self._inhibitors = inhibitors.copy()
    inhibitors = property(_get_inhibitors, _set_inhibitors)

    def _get_data(self):
        return self._measurement
    def _set_data(self, data):
        assert isinstance(data, (int, float))
        self._measurement = data
    data = property(_get_data, _set_data)

    def cues_as_dict(self):
        data = self.stimuli.copy()
        data.update(self.inhibitors)
        return data

    def get_cues(self):
        cues = sorted(self.stimuli.keys()) + sorted(self.inhibitors.keys())
        return cues

    def __str__(self):
        txt = "Protein {} (cellLine {}) measured at time {} has value: {}".format(
                self.protein_name, self.cellLine, self.time, self.data)
        txt += " for the following measurement: \n\tStimuli: {}\n\tInhibitors: {}.".format(
                self.stimuli, self.inhibitors)
        return txt


class Measurements(list):
    """Data structure to store list of measurements

        >>> es = Measurements()
        >>> e1 = Measurement("AKT", 0, {"EGFR":1}, {"AKT":0}, 0.1)
        >>> e2 = Measurement("AKT", 5, {"EGFR":1}, {"AKT":0}, 0.5)
        >>> es.add_measurements([e1,e2])

    """
    def __init__(self, measurements=None):
        if measurements:
            self.add_measurements(measurements)

    def add_measurements(self, measurements):
        if isinstance(measurements, list) is False:
            measurements = [measurements]
        for exp in measurements:
            # FIXME do we need a costly copy here ?
            self.append(copy.deepcopy(exp))

    def _get_species(self):
        species = sorted(list(set([e.protein_name for e in self])))
        return species
    species = property(_get_species)

    def _get_stimuli(self):
        stimuli = set()
        for e in self:
            stimuli = stimuli.union(e.stimuli.keys())
        return stimuli
    stimuli = property(_get_stimuli)

    def _get_inhibitors(self):
        inhibitors = set()
        for e in self:
            inhibitors = inhibitors.union(e.inhibitors.keys())
        return inhibitors
    inhibitors = property(_get_inhibitors)

    def get_protein(self):
        return [x.protein_name for x in self]
    def get_data(self):
        return [x.data for x in self]
    def get_time(self):
        return [x.time for x in self]
    def get_cell(self):
        return [x.cellLine for x in self]



class MIDASBuilder(object):
    """STarts a MIDAS file from scratch and export 2 CSV MIDAS file.

    .. warning:: to be used with care. Right now it seems to work but still in
        development.

        >>> m = MIDASBuilder()
        >>> e1 = Measurement("AKT", 0, {"EGFR":1}, {"AKT":0}, 0.1)
        >>> e2 = Measurement("AKT", 5, {"EGFR":1}, {"AKT":0}, 0.5)
        >>> e3 = Measurement("AKT", 10, {"EGFR":1}, {"AKT":0}, 0.9)
        >>> e4 = Measurement("AKT", 0, {"EGFR":0}, {"AKT":0}, 0.1)
        >>> e5 = Measurement("AKT", 5, {"EGFR":0}, {"AKT":0}, 0.1)
        >>> e6 = Measurement("AKT", 10, {"EGFR":0}, {"AKT":0}, 0.1)
        >>> for e in [e1, e2, e3, e4, e5, e6]:
        ...     m.add_measurement(e)
        >>> m.to_midas("test.csv")

    This class allows one to add measurements to obtain a dataframe compatible with
    XMIDAS class, which can then be saved using XMIDAS.to_midas.

    If an inhibitor or stimuli is not provided, we assume ti is absent (set 
    to zero).


    """
    def __init__(self):
        self.measurements = Measurements()

    def test_example(self, Nspecies=20, N=10, times=[0,10,20,30]):
        """

        N number of stimuli and inhibitors
        Ntime =
        There are duplicates so Nrows = N*2 * Ntimes * 2

        """
        self.measurements = Measurements()
        """e1 = Measurement("DIG1", 0, {"EGF":1, "Akt":1}, {}, 1)
        e2 = Measurement("DIG1", 0, {"EGF":1, "Akt":1}, {}, 1.2)
        e3 = Measurement("DIG1", 0, {"EGF":1, "Akt":1}, {}, 1.3)
        e4 = Measurement("DIG2", 0, {"EGF":1, "Akt":1}, {}, 3)
        e5 = Measurement("DIG2", 10, {"EGF":1, "Akt":1}, {}, 2)
        e6 = Measurement("DIG1", 10, {"EGF":1, "Akt":1}, {}, 1.5)
        e7 = Measurement("DIG3", 20, {"EGF":1, "Akt":0}, {}, 3.5)
        self.add_list_measurements([e1,e2,e3,e4,e5,e6, e7])"""

        species = ['AKT' + str(i) for i in range(1,Nspecies)]

        stimuli = ['S'+str(i) for i in range(1, N)]
        inhibitors = ['I'+str(i) for i in range(1, N)]

        for this in species:
            N1 = random.randint(1, N)
            N2 = random.randint(1, N)
            random.shuffle(stimuli)
            random.shuffle(inhibitors)
            d_sti = dict([(x,1) for x in stimuli[0:N1]])
            d_sti.update(dict([(x,0) for x in stimuli[N1:]]))
            d_inh = dict([(x,1) for x in inhibitors[0:N2]])
            d_inh.update(dict([(x,0) for x in inhibitors[N2:]]))
            for time in times:
                e = Measurement(this, time,
                        d_sti, d_inh, random.random())
                self.add_measurements(e)
                # add duplicated values
                e = Measurement(this, time,
                        d_sti, d_inh, random.random())
                self.add_measurements(e)

    def __len__(self):
        return len(self.measurements)

    def add_measurements(self, measurements):
        self.measurements.add_measurements(measurements)

    def get_colnames(self):
        return self.measurements[0].get_cues()

    def _get_stimuli(self):
        return self.measurements.stimuli
    stimuli = property(_get_stimuli)

    def _get_inhibitors(self):
        return self.measurements.inhibitors
    inhibitors = property(_get_inhibitors)

    def get_df_exps(self):
        stimuli = list(self.stimuli)
        inhibitors = list(self.inhibitors)
        Ns = len(stimuli)
        Ni = len(inhibitors)
        Nrows = len(self)
        N = Ns+Ni

        df = pd.DataFrame(np.zeros(N*Nrows).reshape(Nrows, N), 
                index=range(0,Nrows),
                columns=[['Stimuli']*Ns + ['Inhibitors']*Ni, stimuli + inhibitors])
        df.sortlevel(axis=1, inplace=True)

        # this is the slowest part in the 2 next loops.
        for stimulus in stimuli:
             df.loc[:,('Stimuli', stimulus)] = [x.stimuli.get(stimulus, 0) 
                     for x in  self.measurements]
        for inhibitor in inhibitors:
             df.loc[:,('Inhibitors', inhibitor)] = [x.inhibitors.get(inhibitor, 0) 
                     for x in  self.measurements]

        df['time'] = self.measurements.get_time()
        df['cell'] = self.measurements.get_cell()
        df.reset_index(inplace=True)
        df.rename(columns={'index':'experiment'}, inplace=True)

        # set indexes now based on cell name, time and index, which will need
        # to be renamed as condition
        df.set_index(['cell', 'experiment', 'time'], inplace=True)

        groups = df.groupby(by=list(df.columns.values)).groups

        df = df.drop_duplicates()

        experiment_names = ['experiment_%s' % i for i in 
                range(0, len(groups.keys()))]

        df.reset_index(inplace=True)
        self._df = df.copy()

        # add dummy inhibitors and stimuli columns to create the
        # strucuture, then remove them
        df['Inhibitors', '__dummy__'] = [1] * df.shape[0]
        df['Stimuli', '__dummy__'] = [1] * df.shape[0]
        df = df[['experiment', 'Inhibitors', 'Stimuli']]
        del df['Inhibitors', '__dummy__']
        del df['Stimuli', '__dummy__']

        # drop time and cell info
        # FIXMEL warning raised from call here below
        df.loc[:,'experiment'] = experiment_names
        #df.set_index(['cell', 'experiment', 'time'], inplace=True)
        df.set_index(['experiment'], inplace=True)
        return df, groups

    def get_df_data(self):
        df = pd.DataFrame({
            'protein': self.measurements.get_protein(),
            'data':  self.measurements.get_data(),
            'time': self.measurements.get_time(),
            'cell': self.measurements.get_cell()
            })
        df.reset_index(inplace=True)
        df.rename(columns={'index':'experiment'}, inplace=True)
        return df

    def _get_xmidas(self):
        # FIXME if no time zero provided, assumes this is a fold change 
        # and set all values to zero.

        from cno.io.midas import XMIDAS
        xm = XMIDAS()
        if len(self.measurements) == 0:
            return xm

        xm.df = self.get_df_data()
        df_exps, groups = self.get_df_exps()
        # set the name of the experiments in the first df (df_exps)
        xm._experiments = df_exps

        # set name of the experiments in second df (based on the group)
        # FIXME: multi cell line fails here
        mapping = {}
        for i, name_exp in enumerate(df_exps.index): # loop over experiments
            # looking for the indices (row index) of the experiments
            # within a group
            exps = groups[tuple(df_exps.ix[i].values)]
            exp_index = [exp[1] for exp in exps]
            # all those rows should be named with same experiment name
            # that is the variabe we are looping on (name)
            for j in exp_index:
                mapping[j] = name_exp

            # we can now rename the experiment in the dataframe that containts the data
            # indeed for now experiments are just values from 0 to N
            xm.df.loc[exp_index, 'experiment'] = name_exp

        self._df_data = xm.df.copy()
        # now that experiments have been renamed, we can pivot the dataframe
        xm.df  = pd.pivot_table(xm.df,
                            index=['cell', 'experiment', 'time'],
                            columns='protein', values='data')
        #xm.df.reset_index(inplace=True)
        #xm.df['experiment'] = [mapping[x] for x in xm.df['experiment'].values]
        #xm.df.set_index(['cell', 'experiment', 'time'], inplace=True)

        xm.create_empty_simulation()
        # cell line must be set to one of the cell lines
        xm.cellLine = xm.cellLines[0]

        xm.errors = xm.df * 0

        return xm
    xmidas = property(_get_xmidas)

