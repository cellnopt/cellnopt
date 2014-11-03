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

import numpy as np
import pandas as pd

try:
    from cno.io.midas import XMIDAS
except:
    pass




__all__ = ["Measurement", "Measurements", "MIDASBuilder"]


class Measurement(object):
    """Data structure to store a measurement.

    Givem a list of stimuli and inhibitor, stores a measure
    at a given time.


        >>> from cno.io.midas_extra import Measurement
        >>> m = Measurement("AKT", 0, {"EGFR":1}, {"AKT":0}, 0.1)

    """
    def __init__(self, name, time, stimuli, inhibitors, measurement,
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

        self._measurement = measurement
        self._name = name

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

    def _get_name(self):
        return self._name
    def _set_name(self, name):
        assert isinstance(name, str)
        self._name = name
    name = property(_get_name, _set_name)

    def _get_stimuli(self):
        return self._stimuli
    def _set_stimuli(self, stimuli):
        isinstance(stimuli, dict)
        for k,v in stimuli.iteritems():
            if v >1 or v<0:
                raise ValueError("Value of the stimulus {} must be inside the range [0,1]".format(k))
        self._stimuli = stimuli.copy()
    stimuli = property(_get_stimuli, _set_stimuli)

    def _get_inhibitors(self):
        return self._inhibitors
    def _set_inhibitors(self, inhibitors):
        isinstance(inhibitors, dict)
        for k,v in inhibitors.iteritems():
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
                self.name, self.cellLine, self.time, self.data)
        txt += " for the following measurement: \n\tStimuli: {}\n\tInhibitors: {}.".format(
                self.stimuli, self.inhibitors)
        return txt


class Measurements(object):
    """Data structure to store list of measurements

        >>> es = Measurements()
        >>> e1 = Measurement("AKT", 0, {"EGFR":1}, {"AKT":0}, 0.1)
        >>> e2 = Measurement("AKT", 5, {"EGFR":1}, {"AKT":0}, 0.5)
        >>> es.add_single_measurements([e1,e2])

    """
    def __init__(self):
        self.measurements = []

    def add_single_measurements(self, measurements):
        for exp in measurements:
            import copy
            self.measurements.append(copy.deepcopy(exp))

    def _get_species(self):
        species = sorted(list(set([e.name for e in self.measurements])))
        return species
    species = property(_get_species)

    def __len__(self):
        return len(self.measurements)


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

    More sophisticated builders can be added.


    """
    def __init__(self):
        self.measurements = []

    def test_example(self, Nspecies=20, N=10, times=[0,10,20,30]):
        """

        N number of stimuli and inhibitors
        Ntime =
        There are duplicates so Nrows = N*2 * Ntimes * 2

        """
        self.measurements = []
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

        import random
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
                self.add_measurement(e)
                # add duplicated values
                e = Measurement(this, time,
                        d_sti, d_inh, random.random())
                self.add_measurement(e)

    def __len__(self):
        return len(self.measurements)

    def add_measurement(self, measurement):
        self.measurements.append(measurement)

    def add_list_measurements(self, measurements):
        for e in measurements:
            self.add_measurement(e)

    def get_colnames(self):
        return self.measurements[0].get_cues()

    def _get_stimuli(self):
        stimuli = set()
        for e in self.measurements:
            stimuli = stimuli.union(e.stimuli.keys())
        return stimuli
    stimuli = property(_get_stimuli)

    def _get_inhibitors(self):
        inhibitors = set()
        for e in self.measurements:
            inhibitors = inhibitors.union(e.inhibitors.keys())
        return inhibitors
    inhibitors = property(_get_inhibitors)

    def _get_measurement_name(self, e):
        data = e.stimuli.copy()
        data.update(e.inhibitors)

        # let us build a dataframe corresponding to the measurement. This is a 1-row DF
        mydf = pd.DataFrame(data, index=[0], columns=self._dfexp.columns)

        # let us compare it with the full list of unique measurements to figure out
        # its name
        candidate = self._dfexp[(self._dfexp == mydf.ix[0]).all(axis=1)]
        candidate = candidate.index # just need the indices, not the content
        if len(candidate) != 1:
            print(candidate)
            raise ValueError("Found 0 or more than 1 candidate measurement. ")
        return candidate[0]

    def get_df_exps(self):
        #m.test_example(times=range(0,60))
        # %timeit -n 3 df,dfd = m.get_df_exps()
        # 60 ms per loop

        stimuli = list(self.stimuli)
        inhibitors = list(self.inhibitors)
        Ns = len(stimuli)
        Ni = len(inhibitors)
        Nrows = len(self)
        N = Ns+Ni

        df = pd.DataFrame(np.arange(N*Nrows).reshape(Nrows, N), index=range(0,Nrows),
                columns=[['Stimuli']*Ns + ['Inhibitors']*Ni, stimuli + inhibitors])
        df.sortlevel(axis=1, inplace=True)


        # FIXME: this is the slowest part in the 2 next loops.
        for stimulus in stimuli:
             df.loc[:,('Stimuli', stimulus)] = [x.stimuli[stimulus] for x in
                 self.measurements]
        for inhibitor in inhibitors:
             df.loc[:,('Inhibitors', inhibitor)] = [x.inhibitors[inhibitor] for x in
                 self.measurements]

        df['time'] = [x.time for x in self.measurements]
        df['cell'] = [x.cellLine for x in self.measurements]
        df.reset_index(inplace=True)
        df.rename(columns={'index':'experiment'}, inplace=True)

        # set indexes now based on cell name, time and index, which will need
        # to be renamed as condition
        df.set_index(['cell', 'experiment', 'time'], inplace=True)

        groups = df.groupby(by=list(df.columns.values)).groups

        df = df.drop_duplicates()
        return df, groups

    def get_df_data(self):

        df = pd.DataFrame({
            'protein': [x.name for x in self.measurements],
            'data':  [x.data for x in self.measurements],
            'time': [x.time for x in self.measurements],
            'cell':   [x.cellLine for x in self.measurements]
            })
        df.reset_index(inplace=True)
        df.rename(columns={'index':'experiment'}, inplace=True)

        # NOTE the pivot method does not allow multi index in the pandas
        # version used during the development of this method.
        df = pd.pivot_table(df,
                       index=['cell', 'experiment', 'time'],
                        columns='protein', values='data')

        return df

    def _get_xmidas(self):

        #df_data = self.get_df_data()

        #species = sorted(list(set([e.name for e in self.measurements])))
        #times = sorted(list(set([e.time for e in self.measurements])))
        #measurement_names = list(df_measurements.index)

        #df = self._get_xmidas_df()
        #df = df.sort_index(axis=1)
        #df = df.sortlevel(["experiment"])

        x = XMIDAS()

        # proper column compatible with MIDAS:
        #stimuli = sorted(list(set([y for e in self.measurements for y in e.stimuli.keys()])))
        #inhibitors = sorted(list(set([y for e in self.measurements for y in e.inhibitors.keys()])))
        x.df = self.get_df_data()
        x.create_empty_simulation()
        x._experiments = self.get_df_exps()


        x._cellLine = x.cellLines[0]

        # FIXME do we need this raw attribute
        #x._rawdf = x.copy()
        #x._rawexp = x._measurements.copy()

        return x
    xmidas = property(_get_xmidas)

    def to_midas(self, filename):
        xmidas = self.xmidas
        xmidas.save2midas(filename)







