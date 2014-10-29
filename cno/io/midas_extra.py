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
import types

import pylab
import numpy
import numpy as np
import pandas as pd

try:
    from cno.io.midas import XMIDAS
except:
    pass


from easydev.logging_tools import Logging
import colormap
from easydev import check_param_in_list


__all__ = ["Measurement", "Measurements", "MIDASBuilder"]


class Measurement(object):
    """Data structure to store a measurement.

    Givem a list of stimuli and inhibitor, stores a measure
    at a given time.


        >>> from cno.io.midas_extra import Measurement
        >>> m = Measurement("AKT", 0, {"EGFR":1}, {"AKT":0}, 0.1)

    """
    def __init__(self, protein_name, time, stimuli, inhibitors, measurement,
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
                self.protein_name, self.cellLine, self.time, self.data)
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
        species = sorted(list(set([e.protein_name for e in self.measurements])))
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

    def test_example(self):
        self.measurements = []
        """e1 = Measurement("DIG1", 0, {"EGF":1, "Akt":1}, {}, 1)
        e2 = Measurement("DIG1", 0, {"EGF":1, "Akt":1}, {}, 1.2)
        e3 = Measurement("DIG1", 0, {"EGF":1, "Akt":1}, {}, 1.3)
        e4 = Measurement("DIG2", 0, {"EGF":1, "Akt":1}, {}, 3)
        e5 = Measurement("DIG2", 10, {"EGF":1, "Akt":1}, {}, 2)
        e6 = Measurement("DIG1", 10, {"EGF":1, "Akt":1}, {}, 1.5)
        e7 = Measurement("DIG3", 20, {"EGF":1, "Akt":0}, {}, 3.5)
        self.add_list_measurements([e1,e2,e3,e4,e5,e6, e7])"""

        species = ['AKT' + str(i) for i in range(1,20)]

        N = 10
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
            for time in [0,10,20,30]:

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
        import pandas as pd
        import numpy as np

        stimuli = list(self.stimuli)
        inhibitors = list(self.inhibitors)
        Ns = len(stimuli)
        Ni = len(inhibitors)
        Nrows = len(self)
        N = Ns+Ni

        df = pd.DataFrame(np.arange(N*Nrows).reshape(Nrows, N), index=range(0,Nrows), 
                columns=[['Stimuli']*Ns + ['Inhibitors']*Ni, stimuli + inhibitors])

        df.sortlevel(axis=1, inplace=True)

        for stimulus in stimuli:
             df.loc[:,('Stimuli', stimulus)] = [x.stimuli[stimulus] for x in self.measurements]
        for inhibitor in inhibitors:
             df.loc[:,('Inhibitors', inhibitor)] = [x.inhibitors[inhibitor] for x in self.measurements]

        # Create a key (as a string) from all inhibitors/stimuli values
        DATA = [str(df.ix[i].values) for i,x in enumerate(df.index)]

        # get time/cell and index to build a multi index
        df['time'] = [x.time for x in self.measurements]
        df['cell'] = [x.cellLine for x in self.measurements]
        df.reset_index(inplace=True)

        df['summary'] = DATA
        #df.groupby('summary').groups


        # set indexes
        df = df.set_index(['cell', 'index', 'time'])


        g = df.groupby('summary').groups
        Nexp = len(g)



        # NOW, the data. indexes will be the same as above
        proteins = [x.protein_name for x in self.measurements]
        data = [x.data for x in self.measurements]

        # df for the data
        dfd = pd.DataFrame({'protein':proteins, 'data': data, 'time':time, 'cell': cell})
        dfd = dfd.reset_index().pivot('index','protein')

        #df.duplicated(subset=df[['Stimuli', 'Inhibitors']])
        # measurements
        return df, dfd

    """

    In [85]: df = pd.DataFrame(np.arange(3*18).reshape(3,18), index=[['A','A','A'],['exp1', 'exp1', 'exp2'],
    [0,10,0]], columns=[['Stimuli']*9 + ['Inhibitors']*9, list(m.stimuli) + list(m.inhibitors)])

    In [86]: df.index.names = ['a','b','c']

    In [87]: df.sortlevel(axis=0, inplace=True, level=['b','c'])

    In [88]: df.sortlevel(axis=1, inplace=True)


    """

    def _get_xmidas_df(self):
        # this is a df used to stored and extract information from the
        # measurements
        #exp_names = self._get_measurements()
        self._dfexp = self._get_measurements().copy()
        if "TR:" in self._dfexp.columns[0]:
            self._shift=3
            self._dfexp.columns = [col[3:] for col in self._dfexp.columns] # get rid of TR: for query
        else:
            self._shift = 0
        print("building")
        df = pd.DataFrame({
            'time': [this.time for this in self.measurements],
            #'measurement':[self._get_measurement_name(exp_names, this) for this in exps],
            'measurement':[self._get_measurement_name(this) for this in self.measurements],
            'value':[this.data for this in self.measurements],
            'species':[this.protein_name for this in self.measurements],
            'cellLine':[this.cellLine for this in self.measurements]})

        # create multi index data frame
        df.set_index(["cellLine", "experiment", "time"], inplace=True)
        print(2)

        # now we need to move species names as columns and values column as a
        # matrix. There may be NA values.

        # I thought from here, a simple pivot_table call would make the trick of
        # setting species as the column replicates are averaged, which we do not
        # want.

        # amatrix of values. The tricky part is that measuremets/replicaes may
        # be done for a species and not others so there are possibly NA. The
        # matrix that holds the data are a dimension NxM computed here below

        # values will be populated little by little by appending measurements
        # what we know for sure right now are the list of speices and the first
        # experiment index
        Nspecies = len(df.species.unique())
        tuples = []
        df_values = []
        species = sorted(df.species.unique())

        # FIXME: the usage of the query is slow. need to be speed up
        self._df = df
        for index in df.index.unique():
            data = df.xs(index, level=["cellLine", "experiment", "time"])
            # for each combi of cellline/exp/time, what is the max number of
            # replicates over all species
            M = data.groupby("species").species.count().max()
            for x in range(0,M):
                tuples.append(index)

            values = pd.DataFrame([[np.nan] * Nspecies]*M, columns=species)
            for this in species:
                data_species = data.query("species=='{}'".format(this), engine="python")['value']
                if len(data_species):
                    values[this] =  list(data_species) + [np.nan] * (M-len(data_species))
                else:
                    pass
            df_values.append(values)

            # here we built the dt for a single index
        print(3)
        values = pd.concat(df_values, ignore_index=True)
        index = pd.MultiIndex.from_tuples(tuples, names=["cellLine", "experiment", "time"])
        newdf = pd.DataFrame(values.as_matrix(), index=index, columns=species)
        print(4)
        return newdf

    def _get_xmidas(self):
        """pbl: replicates are ignored !!

        .. todo:: get rid of TR: in the measurements df
        """
        df_measurements = self._get_measurements()

        species = sorted(list(set([e.protein_name for e in self.measurements])))
        times = sorted(list(set([e.time for e in self.measurements])))
        measurement_names = list(df_measurements.index)

        df = self._get_xmidas_df()
        df = df.sort_index(axis=1)
        df = df.sortlevel(["experiment"])

        x = XMIDAS()

        # proper column compatible with MIDAS:
        stimuli = sorted(list(set([y for e in self.measurements for y in e.stimuli.keys()])))
        inhibitors = sorted(list(set([y for e in self.measurements for y in e.inhibitors.keys()])))

        columns = []
        for name in df_measurements.columns:
            if name in stimuli:
                if name.startswith("TR:")==False:
                    newname = "TR:" + name
                    columns.append(newname)
                else:
                    columns.append(name)
            elif name in inhibitors:
                if name.startswith("TR:")==False:
                    newname = "TR:" + name + ":i"
                    columns.append(newname)
                else:
                    columns.append(name)
            else:
                raise ValueError()
        df_measurements.columns = columns

        x._measurements = df_measurements.copy()
        x.df = df.copy()
        x.create_empty_simulation()
        x._cellLine = x.cellLines[0]
        x._rawdf = df.copy()
        x._rawexp = x._measurements.copy()
        return x
    xmidas = property(_get_xmidas)

    def _get_measurements(self):
        df_measurements = pd.DataFrame([e.cues_as_dict() for e in
            self.measurements], columns=self.get_colnames())
        df_measurements = df_measurements.drop_duplicates()
        df_measurements.index = range(0, df_measurements.shape[0])
        df_measurements.index = ["measurement_{}".format(this) for this in  df_measurements.index]
        return df_measurements

    def to_midas(self, filename):
        xmidas = self.xmidas
        xmidas.save2midas(filename)







