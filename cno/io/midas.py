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
import numpy as np
import pandas as pd

from easydev.logging_tools import Logging
import colormap

from measurements import Measurement
import midas_normalisation as normalisation

from cno.core import DevTools
from cno.misc import CNOError


__all__ = ["XMIDAS", "Trend", 'MIDASReader']



class MIDAS(object):
    """This is s a description of the MIDAS format (let us call it 1.0)



    column with those name are ignored:

    ['NOINHIB', 'NOCYTO', 'NOLIG', 'NO-CYTO', 'NO-INHIB', 'NO-LIG']

    """
    _valid_codes = {
                   'ID':'identifier',
                   'TR':'stimuli/inhibitors',
                   'DA':'times',
                   'DV':'measurements'}

    _ignore_codes = ['NOINHIB', 'NOCYTO', 'NOLIG', 'NO-CYTO', 'NO-INHIB', 'NO-LIG']

    def __init__(self, filename=None, verbose=False):
        """.. rubric:: Constructor

        :param str filename: the filename (including correct path)

        """
        self._header = []
        self._names = []
        self.verbose = verbose

        #: filename of the original data
        self.filename = filename

        self.cmap_scale = 1 # 1 is same as CellNOptR
        self.fontsize = 16
        self.logging = Logging(verbose)
        self._colormap = colormap.Colormap()

        self._dev = DevTools()


class MIDASReader(MIDAS):
    """A MIDAS reader that converts the experiments and data into dataframes

    Used by XMIDAS to read and convert the data.

    CellLine must be encoded in the header as follows::

        TR:name:CellLine

    If name is not provided, it is replaced internally with "undefined".

    Those columns are ignored in MIDASReader:

        ['NOINHIB', 'NOCYTO', 'NOLIG', 'NO-CYTO', 'NO-INHIB', 'NO-LIG']

    User should not need to use this class. Use :class:`XMIDAS` instead.

    """
    def __init__(self, filename, verbose='ERROR'):
        super(MIDASReader, self).__init__(filename, verbose)
        self._missing_time_zero = True  # FIXME what is this ?

    def read(self):
        # If the provided filename is not defined, nothing to do
        if self.filename == None:
            return
        self.logging.debug("Reading the data")
        # skip spaces after delimiter
        self._data = pd.read_csv(self.filename, skipinitialspace=True, 
                sep=",")
        self.logging.debug(" - Processing cell lines")
        # figure out the cell line names
        self._preprocess_cellines()

        # remove columns that are invalid and check MIDAS validity
        self._midas_validity()
        self.logging.debug(" - Checking format")

        # some cleanup to remove columns that have to be ignored
        labels = ["TR:"+x for x in self._ignore_codes 
                if "TR:"+x in self._data.columns]

        self._data = self._data.drop(labels, axis=1)
        # just clean up the labels (strip)
        self._data.columns = [x.strip() for x in self._data.columns]

        # from the data, build up the experiment and data dataframes
        self.logging.debug(" - Initialising")
        self._init()
        self.logging.debug(" - Data loaded")

    def _preprocess_cellines(self):
        #CellLine are tricky to handle with the MIDAS format because they use the
        #same prefix TR: as the treatments. You must be sure that (1) there are
        # 2 : signs (2) the suffix is called CellLine and (3) the middle name is
        # provided

        # Those could be recognised but only CellLine is correct. So, let us
        # enforce the correct one instead of dealing with all kind of user
        # choices.
        cellLines = [col for col in self._data.columns if "CellLine" in col]
        if len(cellLines) == 0:
            self.logging.warning("Could not find any column with the required keyword 'CellLine'"
            + "capitalisation matters !. Adding a column named TR:undefined:CellLine")
            self._data['TR:undefined:CellLine'] = [1] * self._data.shape[0]

        for name in cellLines:
            if name.count(":") != 2:
                txt = "Column name related to CellLine %s must have 2 : characters".format(name)
                txt += "{} has less or more than 2"
                raise CNOError(txt)

        # rename the cellLines if there are undefined
        # let us give thme the same name ('undefined')
        # this should be done at the level of the data
        count = 0
        for i, name in enumerate(cellLines):
            if name.split(":")[1] == "":
                self.logging.warning("Found a CellLine column without name. Renamed to undefined.")
                columns = list(self._data.columns)
                columns[i] = "TR:undefined:CellLine"
                self._data.columns = columns
                count += 1
        if count > 1:
            raise CNOError("More than 1 undefined cell line. " +
            "Please name the cellLine columns with TR:name:CellLine")
        # if there is only one undefined, no need for an error.
        # otherwise they will be an ambiguity.

        celltype_names = self._get_cellLines()

        if len(cellLines) >= 2:
            if self.cellLine == None or self.cellLine not in celltype_names:
                txt = "More than 1 celline was found.\n"
                txt += "You must select one amongst: {}".format([this.split(":")[1] for this in cellLines])
                raise CNOError(txt)
            else:
                # we could remove columns and rows where cell type is not correct.
                # but we are going to do it in the _init so that a user can
                # change his mind and reset the cell line to be looked at.
                pass

        if len(celltype_names) == 1:
            self._cellLine = celltype_names[0]

    def _midas_validity(self):

        # check validity of TR:
        for this in self._data.columns:
            if ":" not in this:
                txt = "Error in header of the input MIDAS file\n"
                txt += "    Column's name must contain the special character ':'\n"
                txt += "    Found {0}".format(this)
                raise CNOError(txt)
            if this.split(":")[0] not in self._valid_codes.keys():
                txt = "Error in header of the input MIDAS file\n"
                txt += "    Column's name must start with one of the valid code {}\n".format(self._valid_codes)
                raise CNOError(txt)

        # check if zero time data is available otherwise need to call
        # maybe better to play with -df
        #_duplicate_time_zero_using_inhibitors_only
        df = self._data[[this for this in self._data.columns if this.startswith("DA")]]
        unique_times = list(set(df.as_matrix().flatten()))
        unique_times.sort()
        if len(unique_times) <2:
            raise CNOError("Must contains at least 2 time points including time zero.")
        if 0 not in unique_times or len(unique_times)<=1:
            raise CNOError("You must have zero times in the MIDAS file, that was not found.")
        times = list(df.as_matrix().flatten())
        counter = {}
        for time in unique_times:
            counter[time] = times.count(time)

        self._missing_time_zero = False
        if len(set(counter.keys())) >= 2:
            if counter[0] != counter[unique_times[1]]:
                # call correct function
                self._missing_time_zero = True
            else:
                pass

        if len([x for x in self._data.columns if x.startswith("DV")]) == 0:
            raise CNOError("Header of MIDAS file has no columns starting with DV. expects at least one")


    def _init(self):

        # select only data that matches the cell line choice made by the user.
        cellLine = 'TR:%s:CellLine' % self.cellLine

        # select data for the cell line provided by the user.
        _data = self._data[self._data[cellLine] == 1]

        # and remove all column with the CellLine keyword now
        _data = _data[ [col for col in _data.columns if "CellLine" not in col]]

        #drop ID columns if any
        _data = _data[ [col for col in _data.columns if col.startswith("ID")==False]]

        # figure out the treatments, times and measurements
        df_tr = _data[[this for this in _data.columns if this.startswith("TR")]]
        df_da = _data[[this for this in _data.columns if this.startswith("DA")]]
        df_dv = _data[[this for this in _data.columns if this.startswith("DV")]]

        # let us gather experiments nd replace empty fields with zeros ??
        # FIXME filling NA with zero is it correct ?
        value_experiments = _data[df_tr.columns].copy()
        value_experiments.replace("NaN", 0, inplace=True)

        value_signals = _data[df_dv.columns].as_matrix()
        value_times = _data[df_da.columns]

        names = [this for this in df_tr[:] if "CellLine" not in this]
        self._experiments = _data[names].drop_duplicates()
        self._experiments.index = range(0, self._experiments.shape[0])
        self._experiments.index = ["experiment_{}".format(this)
            for this in  self._experiments.index]

        # FIXME: already done above ...
        self._experiments.replace("NaN", 0, inplace=True)

        # build the tuples that will be used by the MultiIndex dataframe
        tuples = []
        # make sure to read each row in the original data using .shape
        for ix in _data.index:
            # 1. find name of the experiment by
            this_exp = value_experiments.ix[ix]
            # scan unique experiments and figure out which one is this_exp
            exp_name = None
            for this_unique_exp in self._experiments.iterrows():
                if all(this_unique_exp[1] == this_exp):
                    #[1:] to ignore the cell line TODO
                    # found it
                    exp_name = this_unique_exp[0]
            assert exp_name != None

            # 2. times
            time = set(value_times.ix[ix])
            assert len(time) == 1
            time = list(time)[0]
            tuples.append((self.cellLine, exp_name, time))

        # replace empty strings with 0
        # note that ,,, is interpreted as ,NaN,NaN,NaN
        # but , , , is interpreted as ," "," "," ",
        self._experiments = self._experiments.applymap(lambda x: 0
            if isinstance(x, basestring) and x.isspace() else x)

        # must convert read data into numeric value.
        self._experiments = self._experiments.convert_objects(convert_numeric=True,
                                                              copy=True)

        index =  pd.MultiIndex.from_tuples(tuples, names=self._levels)
        #keep = [this for this in self.df.columns if this not in ["experiments", "times"]]
        names_species = [x for x in _data.columns if x.startswith('DV:')]
        names_species = [x[3:] for x in names_species]

        self.df = pd.DataFrame(value_signals, index=index, columns=names_species)
        self.df = self.df.sortlevel(["experiment"])
        self.df = self.df.sort_index(axis=1) # sort the species

        # Get rid of TR in experiments
        self._experiments.columns = [this.replace("TR:", "") for this in self._experiments.columns]


        cues = []
        for c in self._experiments.columns:
            if c.count(":") >=2:
                raise CNOError("Invalid header. Found more than 2 : sign in a column")

            if c.endswith(":Stimuli"):
                cues.append(c.split(":")[0])
            elif c.endswith(":Inhibitors"):
                cues.append(c.split(":")[0] + ":i")
            elif c.endswith("i"):
                cues.append(c[0:-1] + ":i")
            else:
                cues.append(c)
        self._experiments.columns = cues

        inh = [x for x in self._experiments.columns if x.endswith(":i") is True]
        stim = [x for x in self._experiments.columns if x.endswith(":i") is False]

        # Finally, build up the dataframe for experiments
        self._experiments = pd.DataFrame(self._experiments.values,
                     columns=[['Stimuli']*len(stim) + ['Inhibitors']*len(inh),
                              [x.replace(":i","") for x in self._experiments.columns]],
                        index= self._experiments.index)

        self._experiments.sortlevel(axis=1, inplace=True)

        if self._missing_time_zero == True:
            try:
                self._duplicate_time_zero_using_inhibitors_only()
            except Exception as e:
                print(e.message)
                self.logging.warning("Issue with missing time zero ")
                return

        if self.df.shape[0] > len(self.times) * self.experiments.shape[0]:
            self.logging.warning("you may have duplicated experiments" +
            "please average the replicates using self.average_replicates(inplace=True)")

        if self.df.max(skipna=True).max(skipna=True) > 1:
            self.logging.warning("values larger than 1. " +
            "You may want to normalise/scale the data")


    def _duplicate_time_zero_using_inhibitors_only(self):
        """
        Sometimes the time zero data sets are not explicitly written in MIDAS
        files. One example is MD-ExtLiverHepG2-MCP2010-mod4.csv

        A set of perturbations with same inhibitors but different stimuli
        have the same control, which is the perturbation with same inhibitor
        and no stimuli. That control can therefore be duplicated.

        """
        self.logging.warning("duplicating time zeros data to fit experiment at other times")
        # first figure out the inhibitors

        # find experiment that have missing time zero data
        experiments = list(self.experiments.index)

        tobeadded = None
        for i, this_exp in enumerate(experiments):
            times = self.df.ix[self.cellLine].ix[experiments[i]].index
            if 0 in times:
                pass
            else:
                # need to find the experiment with same inhibitors
                # and time 0
                # there should be only one ? maybe not if replicates.
                these_inhibitors =  self.experiments.ix[this_exp].Inhibitors

                for this_exp_intern in experiments:
                    times = self.df.ix[self.cellLine].ix[this_exp_intern].index
                    if 0 in times:
                        if all(self.experiments.ix[this_exp_intern].Inhibitors ==
                            these_inhibitors):
                            break # so that if there are several replicates found,
                            # we pick up the first one only

                # get the times for this experimenti. it must contain the time zero
                newdata = self.df.xs((self.cellLine, this_exp_intern))
                # we only need the time 0
                newrow = newdata[newdata.index == 0].copy()
                # let us add some index information that is now missing

                newrow[self._levels[0]] = self.cellLine
                newrow[self._levels[1]] = this_exp
                newrow[self._levels[2]] = 0

                # we can now merge with the full data set
                if types.NoneType == type(tobeadded):
                    tobeadded = newrow.copy()
                else:
                    tobeadded = pd.concat([tobeadded, newrow])
        # finally concatenate all the rows with the dataframe df
        df = pd.concat([self.df.reset_index(), tobeadded], ignore_index=True)
        df = df.set_index(self._levels)
        df = df.sortlevel([self._levels[1]]) # experiment
        self.df = df.copy()


class XMIDAS(MIDASReader):
    """The extended MIDAS data structure.

    .. plot::
        :include-source:
        :width: 80%

        from cno import XMIDAS, cnodata
        m = XMIDAS(cnodata("MD-ToyPB.csv"))
        m.df            # access to the data frame
        m.experiments   # access to experiments
        m.sim           # access to simulation
        m.plot()

    Please see main documentation for extended explanation and documentation
    of the methods here below.

    .. warning:: if there are replicates, call average_replicates before creating a
        simulation or calling plot"mode="mse")

    When reading a MIDAS file with several cell lines, all data is stored but
    only one cell line can be visualised and manipulated at a time::

        # when reading, you must give the cell line name you want to activate
        multiple_cell = XMIDAS('somedata', 'cell1')
        # valid cell line names are stored in this attribute
        multiple_cell.cellLines
        # activating another cell line needs to set this attribute:
        multiple_cell.cellLine = 'cell2'

    """
    def __init__(self, filename=None, cellLine=None, verbose=False):
        """**Constructor**

        :param str filename: filename of a MIDAS file or a XMIDAS instance
        :param str cellLine: name of a cell Line (compulsary if several
            cell lines are present)

        """
        super(XMIDAS, self).__init__(filename, verbose=verbose)

        self._cellLine = cellLine

        # multi index related (position of the columns)
        self._celltype_index = 0
        self._experiment_index = 1
        self._time_index = 2
        # names of the multi index level (rows)
        self._levels = ["cell", "experiment", "time"]

        self.verbose = verbose

        #self._ignore_invalid_columns = True

        self._data = pd.DataFrame()
        self._experiments = pd.DataFrame()
        self.df = pd.DataFrame()
        self.sim = self.df.copy() * 0
        self.errors = self.df.copy() * 0

        if filename is None:
            pass
        elif isinstance(filename, str):
            self.read()
            self.create_empty_simulation()
            self.errors = self.sim.copy()
        else:
            try:
                md = filename
                self.filename = md.filename
                self.df = md.df.copy()
                self._ = md._data.copy()
                self._missing_time_zero = md._missing_time_zero
                self.df = self.df.copy()
                self._experiments = md.experiments.copy()
                self._cellLine = md.cellLine
                self.sim = md.sim.copy()
                self.errors = md.errors.copy()
            except Exception as e:
                raise CNOError(e.message +
                "Invalid input file. Expecting a XMIDAS instance or valid filename")
        self._params = {
                'plot.layout.space': 0.5,
                'plot.fontsize':16,
                'plot.fontsize.times':16,
                'plot.fontsize.species':16,
                'plot.fontsize.perturbations':14,
                'plot.fontsize.titles':16,
                'plot.layout.shift_top': 1,
                'plot.colorbar.N': 10,
                }

    def _manage_replicates(self):
        """

        """
        groups = self.df.groupby(level=self._levels).groups
        if any([len(this)>1 for this in groups.values()])==False:
            self.logging.info("No replicates found")
        else:
            newdf = self.df.groupby(level=self._levels).agg([np.mean, np.std])
            return newdf

    def average_replicates(self, inplace=False):
        """Average replicates if any

        :param bool inplace: default to False

        If inplace, a new dataframe :attr:`errors` is created
        and contains the errors (standard deviation)


        """
        df = self._manage_replicates()

        if isinstance(df, pd.DataFrame):
            dfstd = df[[this for this in df.columns if "std" in this]]

            df = df[[this for this in df.columns if "mean" in this]]
            df.columns = [this[0] for this in df.columns]
            if inplace:
                self.df = df.copy()
                self.sim = self.df.copy()
                self.errors = dfstd.copy()
            else:
                return df

    def __div__(self, this):
        m = self.copy()
        m.df /= this
        return m

    def __mul__(self, this):
        m = self.copy()
        m.df *= this
        return m

    def __sub__(self, this):
        # m1 + m2   or m1+=m2
        m = self.copy()
        if isinstance(this, XMIDAS):
            m.df -= this.df
        else:
            m.df -= this
        return m

    def __add__(self, this):
        # m1 + m2   or m1+=m2
        m = self.copy()
        if isinstance(this, XMIDAS):
            m.df += this.df
        else:
            m.df += this
        return m

    def __getitem__(self, item):
        """Simple alias to the dataframe's columns

        .. doctest::

            >>> from cno import XMIDAS, cnodata
            >>> m = XMIDAS(cnodata("MD-ToyPB.csv"))
            >>> m['p38'].ix[0]
            0.1076
            >>> # equivalent to m.df['p38'].ix[0]

        .. note:: this is a read-only access
        """
        if len(item)==3 and isinstance(item, str) != True:
            i1, i2, i3 = item
            return self.df.ix[i1].ix[i2].ix[i3]
        elif isinstance(item,str)==True:
            if item in self.df.columns:
                return self.df[item].copy()

    def _get_experiments(self):
        return self._experiments
    experiments = property(_get_experiments,
            doc="Return dataframe with experiments")

    def _get_names_species(self):
        return self.df.columns.values
    names_species = property(_get_names_species,
            doc="list of species")
    names_signals = property(_get_names_species,
            doc="same as :attr:`names_species`")
    species = property(_get_names_species,
            doc="Getter for the columns of the dataframe that represents the species/signals")
    signals = property(_get_names_species,
            doc="Getter for the columns of the dataframe that represents the species/signals")

    def _get_names_cues(self):
        cues = self.names_stimuli + self.names_inhibitors
        return cues
    names_cues = property(_get_names_cues,
            doc="Return list of stimuli and inhibitors together")

    def _get_cellLines(self):
        # XMIDAS created from an input file
        if len(self._data):
            names = [this for this in self._data.columns if "CellLine" in this]
            names = [this.split(":")[1] for this in names]
            return names
        # XMIDAS created from another builder e.g. MIDASBuilder
        else:
            names = self.df.index.levels[0]
            return names
    cellLines = property(_get_cellLines, doc="Return available cell lines")

    def _get_cellLine(self):
        return self._cellLine
    def _set_cellLine(self, name):
        if name not in self.cellLines:
            raise CNOError("Invalid cellLine name {}. Valid ones are {}".format(name,
                           self.cellLines))
        self._cellLine = name
        # TODO: do we need to call _init again ?
        self._init()
    cellLine = property(_get_cellLine, _set_cellLine,
                        doc="Getter/Setter of the active cell line")

    def _get_times(self):
        times = self.df.index.levels[self._time_index]
        return sorted(list(times))
    times = property(_get_times, doc="Getter to the different times")

    def _get_names_inhibitors(self):
        try:
            return list(self.experiments.Inhibitors.columns)
        except:
            return []
    names_inhibitors = property(_get_names_inhibitors,
                                doc="return list of inhibitors")

    def _get_names_stimuli(self):
        try:
            return list(self.experiments.Stimuli.columns)
        except:
            return []
    names_stimuli = property(_get_names_stimuli,
                             doc="returns list of stimuli")

    def reset(self):
        """Reset the dataframe to original data set"""
        self._init()

    def _check_consistency_data(self):
        # times consistency all times must have same length
        #
        pass

    def reset_index(self):
        """Remove all indices (cellLine, time, experiment)

        Done in the 3 dataframes :attr:`df`, :attr:`sim` and :attr:`errors`
        """
        self.df.reset_index(inplace=True)
        self.sim.reset_index(inplace=True)
        self.errors.reset_index(inplace=True)

    def set_index(self):
        """Reset all indices (cellLine, time, experiment)

        Done in the 3 dataframes :attr:`df`, :attr:`sim` and :attr:`errors`
        """
        if len(self.df):
            self.df.set_index(self._levels, inplace=True)
            self.sim.set_index(self._levels, inplace=True)
            self.errors.set_index(self._levels, inplace=True)

    def remove_species(self, labels):
        """Remove a set of species

        :param labels: list of Species (list of strings) or just one
            species(single string or as a list.

        ::

            m.remove_species("p38")
            m.remove_species(["p38"])

        """
        labels = self._dev.tolist(labels)
        columns = self.df.columns[:]
        for label in labels:
            if label not in columns:
                self.logging.warning("{} not in the species. skipped".format(label))
            else:
                self.df.drop(label, axis=1, inplace=True)
                self.sim.drop(label, axis=1, inplace=True)
                self.errors.drop(label, axis=1, inplace=True)

    def remove_cellLine(self, labels):
        """Remove a cellLine from the dataframe.

        Does not really work since there is only one cellLine in the dataframe.
        If you call this method with the active cell line, the dataframe will
        be empty. You can set it back by activating a new cell line.

        Valid cellLine found in the :attr:`cellLines` attribute.

        This method may be useful later with multi cell line.

        """
        self._remove_labels_from_level(labels, self._levels[0])

    def remove_times(self, labels):
        """Remove time values from the data

        :param list labels: one time point or a list of time points.
            Valid time points are in the :attr:`times` attribute.

        """
        self._remove_labels_from_level(labels, "time")

    def remove_experiments(self, labels):
        """Remove experiment(s) from the dataframe


        :param list labels: one experiment or a list of experiments.
            Valid experiments are in the :attr:`experiments.index` dataframe.
            Experiments are of the form "experiment_12". You can refer
            to an experiment by its number (e.g., here 12).
        """
        self._remove_labels_from_level(labels, "experiment")

    def _remove_labels_from_level(self, labels, level):
        """For a given level and a list of labels, this function
        figures out the rows to remove in the sim/errors/df data frames.

        This is a bit complicated but looks like the proper way

        """
        labels = self._dev.tolist(labels)
        if level == "experiment":
            labels = [x if "experiment" in str(x) else "experiment_"+str(x)
                    for x in labels]

        self.reset_index()
        for label in labels:
            if label not in set(self.df[level]):
                self.logging.warning("{} not found. Skipped".format(label))
            else:
                self.df.drop(self.df[self.df[level] == label].index, inplace=True)
                self.sim.drop(self.sim[self.sim[level] == label].index, inplace=True)
                self.errors.drop(self.errors[self.errors[level] == label].index, inplace=True)
                if level == "experiment":
                    self.experiments.drop(label, inplace=True)
        # and now set the level back as before
        self.set_index()

    def remove_stimuli(self, labels):
        """Remove a stimuli from the :attr:`experiment` dataframe

        :param labels: a string or list of string representing the stimuli
        """
        self._remove_column_experiment(labels, level='Stimuli')

    def remove_inhibitors(self, labels):
        """Remove inhibitor(s) from the :attr:`experiment` dataframe

        :param labels: a string or list of string representing the inhibitor(s)
        """
        self._remove_column_experiment(labels, level='Inhibitors')

    def _remove_column_experiment(self, labels, level):
        if level == 'Inhibitors':
            valid_names = self.names_inhibitors
        elif level == 'Stimuli':
            valid_names = self.names_stimuli
        labels = self._dev.tolist(labels)
        for label in labels:
            # m.experiments.columns.levels[1] contains all names even those that
            # are removed so need to use names_inhibitors instead
            if label not in set(valid_names):
                self.logging.warning("{} not found. Skipped".format(label))
            else:
                self.experiments.drop((level,label), axis=1, inplace=True)

    def rename_stimuli(self, names_dict):
        """Rename stimuli in the :attr:`experiment` dataframe

        :param names_dict: a dictionary with names (keys) to be replaced (values)

        ::

            from cno import XMIDAS, cnodata
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.rename_species({"erk":"ERK", "akt":"AKT"})

        .. seealso:: :meth:`rename_stimuli`, :meth:`rename_inhibitors`

        """
        self._dev.check_param_in_list(names_dict.keys(), self.names_stimuli)
        self.experiments.rename(columns=names_dict, inplace=True)

    def rename_inhibitors(self, names_dict):
        """Rename inhibitors

        :param names_dict: a dictionary with names (keys) to be replaced (values)

        ::

            from cno import XMIDAS, cnodata
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.rename_species({"raf":"RAF"})

        .. seealso:: :meth:`rename_stimuli`, :meth:`rename_species`


        """
        self._dev.check_param_in_list(names_dict.keys(), self.names_inhibitors)
        self.experiments.rename(columns=names_dict, inplace=True)

    def rename_species(self, names_dict):
        """Rename species in the main :attr:`df` dataframe

        :param names_dict: a dictionary with names (keys) to be replaced (values)

        ::


            from cno import cnodata, XMIDAS
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.rename_species({"erk":"ERK", "akt":"AKT"})


        .. seealso:: :meth:`rename_stimuli`, :meth:`rename_inhibitors`
        """
        self._dev.check_param_in_list(names_dict.keys(), self.df.columns)
        columns = list(self.df.columns)
        columns = [c if c not in names_dict.keys() else names_dict[c] for c in columns]
        self.df.columns = columns
        self.sim.columns = columns

    def rename_cellLine(self, to_replace):
        """Rename cellLine indices

        :param dict to_replace: dictionary with mapping of values to be replaced.

        For example; to rename a valid cell line use::

            m.rename_cellLine({"undefined": "PriHu"})

        """
        self.reset_index()
        self.df.replace({self._levels[0]: to_replace}, inplace=True)
        self.set_index()

    def rename_time(self, to_replace):
        """Rename time indices

        :param dict to_replace: dictionary with mapping of values to be replaced.

        For example; to convert time in minutes to time in seconds, use
        something like::

            m.rename_time({0:0,1:1*60,5:5*60})

        """
        self.reset_index()
        # key may be same as value, in which case replace method fails. So, we get rid of them
        to_replace = dict([(k,v) for k,v in to_replace.iteritems() if k!=v])
        self.df.replace({"time": to_replace}, inplace=True)
        self.set_index()

    def merge_times(self, how="mean"):
        """Not implemented yet"""
        raise NotImplementedError

    def add_experiment(self, e):
        """Not implemented yet"""
        raise NotImplementedError

    def scale_max(self, inplace=True):
        r"""Divide all data by the maximum over the entire data set

         .. math::

            X = \frac{X}{M}

        where :math:`M = \max_{e,s,t} X` (with :math:`e` the experiment,
        :math:`s` the species, and :math:`t` the time).
        """
        M = self.df.max().max()
        if inplace:
            self.df /= M
        else:
            return self.df / M

    def scale_max_across_experiments(self, inplace=True):
        r"""Rescale each species column across all experiments


        .. math::

            X_s = \frac{X}{M_s}

        In the MIDAS plot, this is equivalent to dividing each column by
        the max over that column. So, on each column, you should get 1 max values
        set to 1 (if the max is unique). The minimum values may not be set to 0.

        """
        newdf = self.df.divide(self.df.max(), level="experiment")
        if inplace:
             self.df = newdf.copy()
        else:
             return newdf

    def scale_min_max_across_experiments(self, inplace=True):
        r"""Rescale each species column across all experiments

        .. math::

            X_s = \frac{X-m_s}{M-m_s}

         where :math:`m = min_{e,s,t} X` and :math:`M = max_{e,s,t} X`,
        with :math:`e` the experiment, with :math:`s` the species,
        with :math:`t` the time.

        """
        m = self.df.min()
        M = self.df.max()
        data = (self.df - m)/(M-m.astype(np.float64))
        if inplace:
            self.df = data.copy()
        else:
            return data

    def scale_min_max(self, inplace=True):
        r"""Divide all data by the maximum over entire data set

        .. math::

            X = \frac{X-m}{M-m}

        where :math:`m = min_{e,s,t} X` and :math:`M = max_{e,s,t} X`,
        with :math:`e` the experiment, with :math:`s` the species,
        with :math:`t` the time.

        This is an easy (but naive way) to set all data points between 0 and 1.
        """
        m = self.df.min().min()
        M = self.df.max().max()

        if inplace:
            self.df -= m
            self.df /= (M-m)
        else:
            newdf = self.df - m
            newdf /= (M-m)
            return newdf

    def create_empty_simulation(self):
        """Populate the simulation dataframe with zeros.

        The simulation has the same layout as the experiment.

        The dataframe is stored in :attr:`sim`.

        """
        self.sim = self.df * 0

    def create_random_simulation(self):
        """Populate the simulation dataframe with uniformly random values.

        The simulation has the same layout as the experiment.

        The dataframe is stored in :attr:`sim`.
        """
        self.sim = self.df * 0 + np.random.uniform(size=self.df.shape)

    def get_diff(self, sim=None, squared=True, normed=True):
        r"""Return dataframe with differences between the data and simulation

        The dataframe returned contains the MSE (mean square error) by default.

        :param sim: if not provided, uses the :attr:`sim` attribute.
        :param normed:
        :param bool square: set to True to get MSE otherwise, returns the absolute
            values without squaring


        if square is True and normed is True:

        .. math::

            \epsilon = \frac{1}{N_t} \left(X - X^s\right)^2

        We then sum over time and if *normed* is False, :math:`N_t` is set to 1.

        """
        # dataframe cannot be compared to None so, we need this trick:
        if isinstance(sim, types.NoneType):
            sim = self.sim

        if squared is True:
            diff = (sim - self.df).abs()**2
        else:
            diff = (sim - self.df).abs()

        diff = diff.sum(level="experiment")

        if normed:
            N = len(self.times)
            diff = diff/float(N)
        return diff

    # PLOTTING
    def corr(self, names=None, cmap='gist_heat_r', show=True):
        """plot correlation between the measured species

        :param list names: restriction to some species if provided.
        :param string cmap: a valid colormap (e.g. jet). Can also use "green" or "heat".

        .. plot::
            :include-source:
            :width: 80%

            >>> from cno import XMIDAS, cnodata
            >>> m = XMIDAS(cnodata("MD-ToyPB.csv"))
            >>> m.corr(cmap="green")

        """
        cmap = self._get_cmap(cmap)
        corr = self.df.corr()

        if show is True:
            N = corr.shape[0]
            names = self.df.columns[:]
            pylab.clf()
            pylab.pcolor(corr, edgecolors="k", cmap=cmap);
            pylab.xticks([x+0.5 for x in range(0,N)], names, rotation=90)
            pylab.yticks([x+0.5 for x in range(0,N)], names)
            pylab.tight_layout()
            pylab.colorbar()
        return corr

    def plot_sim_data(self, markersize=3, logx=False, linestyle="--", lw=1,
            color="b", marker="x", **kargs):
        """plot experimental curves

        .. plot::
            :width: 80%
            :include-source:

            >>> from cno import XMIDAS, cnodata
            >>> m = XMIDAS(cnodata("MD-ToyPB.csv"));
            >>> m.plot_layout()
            >>> m.plot_data()
            >>> m.plot_sim_data()

        """
        times = np.array(self.times)
        # if simulation do not have the same number of points as data
        simtimes = np.array(self.sim.index.levels[2])

        if logx == False:
            # a tick at x = 0, 0.5 in each box (size of 1) + last x=1 in last box
            xt = pylab.linspace(0, self.nSignals, self.nSignals*2+1)
            M = max(max(times), max(simtimes))
            times = times/float(max(times))
            simtimes = simtimes/float(M)
            xtlabels = self._get_xtlabels()
        else:
            M = float(max(times))
            xtlin = pylab.linspace(0, self.nSignals, self.nSignals*2+1)
            xt = [int(x)+pylab.log10(1+pylab.mod(x,1)*M)/pylab.log10(1+M)
                    for i,x in enumerate(xtlin)]
            xtlabels = self._get_xtlabels()
            M = max( max(pylab.log10(1+times)), max(pylab.log10(1+simtimes)))
            times = pylab.log10(1+times)/max(pylab.log10(1+times))
            simtimes = pylab.log10(1+simtimes)/float(M)
        #for isim, sim in enumerate(self.sim):

        for i in range(0, self.nExps):
            for j in range(0, self.nSignals):
                # divide data by 1.1 to see results close to 1.
                signal = self.names_signals[j]
                exp = self.experiments.index[i]
                data = self.sim[signal][self.cellLine][exp]

                # sometimes we may want to get rid of all NA and show the lines.
                data = data.dropna()
                times = np.array(list(data.index))
                simtimes = times/float(M)

                pylab.plot(simtimes+j, data/1.05+(self.nExps-i-1), marker=marker, color=color,
                        linestyle=linestyle,
                           markersize=markersize, lw=lw)
                #    plot(times+j, sim[i,j]/1.05+(self.nExps-i-1), 'b--o',
                #           markersize=markersize)
        pylab.gca().set_xticklabels(xtlabels, fontsize=kargs.get("fontsize", 10))
        pylab.gca().set_xticks(xt)

    def plot_data(self, markersize=4, logx=False,color="black", **kargs):
        """plot experimental curves

        .. plot::
            :width: 80%
            :include-source:

            >>> from cno import XMIDAS, cnodata
            >>> m = XMIDAS(cnodata("MD-ToyPB.csv"));
            >>> m.plot_layout()
            >>> m.plot_data()

        .. note:: called by :meth:`plot`
        .. seealso:: :meth:`plot`, :meth:`plot_layout`, :meth:`plot_sim_data`
        """
        mode = kargs.get("mode", "data")
        times = np.array(self.times)
        max_time = float(max(self.times))

        if self.df.max().max() > 1 or self.df.min().min() < 0:
            raise ValueError('data must be between 0 and 1')



        if logx == False:
            # a tick at x = 0, 0.5 in each box (size of 1) + last x=1 in last box
            xt = pylab.linspace(0, self.nSignals, self.nSignals*2+1)
            times = times/max_time
            xtlabels = self._get_xtlabels(logx=True)
        else:
            M = max_time
            xtlin = pylab.linspace(0, self.nSignals, self.nSignals*2+1)
            xt = [int(x)+pylab.log10(1+pylab.mod(x,1)*M)/pylab.log10(1+M)
                    for i,x in enumerate(xtlin)]
            xtlabels = self._get_xtlabels()
            times = pylab.log10(1+times)/max(pylab.log10(1+times))

        trend = Trend()

        # TODO must be using the index instead of a range ince indices may not
        # start at zero
        for i in range(0, len(self.experiments)):
            #vMax = float(self.df.max(skipna=True).max(skipna=True))
            for j in range(0, self.nSignals):
                y = self.df[self.names_species[j]][self.cellLine][self.experiments.index[i]]
                
                trend.set(y)
                Y = y

                ratio = 0.95
                if mode == "data":
                    color = trend.get_bestfit_color()

                    pylab.plot(trend.normed_times+j,
                               ratio*Y + self.nExps-i-1 ,
                               'k-o', markersize=markersize, color=color, mfc='gray')

                    pylab.fill_between(trend.normed_times+j,
                                       ratio*Y + self.nExps-1-i ,
                                       self.nExps-1-i, alpha=trend.alpha/1.1,
                                       color=color)
                else:
                    pylab.plot(trend.normed_times+j,
                               ratio*Y + self.nExps-i-1 , 'k-o',
                               markersize=markersize, color="k", mfc='gray')

                #    plot(times+j, sim[i,j]/1.05+(self.nExps-i-1), 'b--o', markersize=markersize)
        pylab.gca().set_xticklabels(xtlabels, fontsize=kargs.get("fontsize",10))
        pylab.gca().set_xticks(xt)

    def _get_cmap(self, cmap=None):
        if cmap == "heat":
            cmap = self._colormap.get_cmap_heat_r()
        elif cmap == "green":
            cmap = self._colormap.get_cmap_red_green()
        return cmap


    def _get_diff_data(self, mode, squared=True):
        diffs = self.get_diff(self.sim, squared=squared)
        # we re-order the dataframe with the indices of the experiments
        # we will do the same with the data !
        diffs = diffs.ix[self.experiments.index]

        if mode == "mse":
            diffs = np.ma.array (diffs, mask=np.isnan(diffs))
        elif mode == "data":
            pass
        return diffs

    def plot_layout(self, cmap="heat",
        rotation=90, margin=0.05, colorbar=True, vmax=None, vmin=0.,
        mode="data", **kargs):
        """plot MSE errors and layout

        :param cmap:
        :param rotation:
        :param margin:
        :param colorbar:
        :param vmax:
        :param vmin:
        :param mode:

        .. plot::
            :width: 80%
            :include-source:

            >>> from cno import cnodata, XMIDAS
            >>> m = XMIDAS(cnodata("MD-ToyPB.csv"));
            >>> m.plot_layout()

        .. note:: called by :meth:`plot`
        .. seealso:: :meth:`plot`, :meth:`plot_layout`, :meth:`plot_sim_data`

        """
        # Get the cmap
        if mode == "data":
            #should be one with zero being white
            cmap = self._colormap.get_cmap_heat_r()
        else:
            cmap = self._get_cmap(cmap)

        # The data to show in the main subplot
        diffs = self._get_diff_data(mode)

        # aliases
        fct = self._params['plot.fontsize.times']
        fcs = self._params['plot.fontsize.species']
        fcp = self._params['plot.fontsize.perturbations']
        fcti = self._params['plot.fontsize.titles']

        # The layout
        gs = pylab.GridSpec(10, 10,
                wspace=self._params['plot.layout.space'],
                hspace=self._params['plot.layout.space'])
        fig = pylab.figure(num=1, figsize=(10, 6))
        shift_top = self._params['plot.layout.shift_top']  #2
        layout_width = 7
        w1 = layout_width + 1
        w2 = layout_width + 2

        ax_main = fig.add_subplot(gs[shift_top:, 0:layout_width])
        ax_main.set_yticks([], [])
        ax_main.set_xticks([], [])

        # stimuli
        ax_stim = fig.add_subplot(gs[shift_top:, layout_width:w1])
        ax_stim.set_yticks([], [])
        ax_stim.set_xticks([], [])
        ax_stim_top = fig.add_subplot(gs[0:shift_top, layout_width:w1])
        ax_stim_top.set_yticks([], [])
        ax_stim_top.set_xticks([], [])

        # inhibitors
        ax_inh = fig.add_subplot(gs[shift_top:, w1:w2])
        ax_inh.set_yticks([], [])
        ax_inh.set_xticks([], [])
        ax_inh_top = fig.add_subplot(gs[0:shift_top, w1:w2])
        ax_inh_top.set_yticks([], [])
        ax_inh_top.set_xticks([], [])

        # colorbar
        if colorbar and mode == 'mse':
            ax_cb = fig.add_subplot(gs[shift_top:, w2:10])
            #ax_cb.set_yticks([], [])
            ax_cb.set_xticks([], [])

        # MAIN subplot with signals
        M = np.nanmax(diffs) # figure out the maximum individual MSE
        m = np.nanmin(diffs) # figure out the minimum individual MSE

        vmax_user = vmax
        vmax = max(1, M)         # if M below 1, set the max to 1 otherwise to M
        if vmax_user:
            vmax = vmax_user

        pylab.sca(ax_main)
        if mode == "mse":
            try: cmap.set_bad("grey", 1.)
            except: pass
            pylab.pcolormesh(pylab.flipud(diffs)**self.cmap_scale, cmap=cmap,
                             vmin=vmin, vmax=vmax, edgecolors='k');
        elif mode == "data":
            try:cmap.set_bad("grey", 1.)
            except: pass
            pylab.pcolormesh(pylab.flipud(diffs*0), cmap=cmap, edgecolors='k');

        # Some tuning of the main plot
        pylab.sca(ax_main)
        pylab.axis([0, diffs.shape[1], 0, diffs.shape[0]])

        # the species name on the top
        ax2 = ax_main.twiny()
        ax2.set_xticks([i+.5 for i,x in enumerate(self.names_species)])
        N = len(self.names_species)
        ax2.set_xticks(pylab.linspace(0.5, N-1, N))
        ax2.set_xticklabels(self.names_species, rotation=90,
                fontsize=fcs)

        # the stimuli
        if len(self.names_stimuli)>0:
            pylab.sca(ax_stim)
            stimuli = np.where(np.isnan(self.stimuli)==False, self.stimuli, 0.5)
            pylab.pcolor(1-pylab.flipud(stimuli), edgecolors='gray', cmap='gray',
                     vmin=0,vmax=1);
            #ax_stim.set_yticks([], [])
            ax_stim.set_xticks([i+.5 for i,x in enumerate(self.names_stimuli)])
            ax_stim.set_xticklabels(self.names_stimuli, rotation=rotation,
                    fontsize=fcp)
            pylab.axis([0,self.stimuli.shape[1], 0, self.stimuli.shape[0]])

        # the inhibitors
        if len(self.names_inhibitors)>0:
            pylab.sca(ax_inh)
            inhibitors = np.where(np.isnan(self.inhibitors)==False,
                                  self.inhibitors, 0.5)
            pylab.pcolor(1-pylab.flipud(inhibitors), edgecolors='gray',
                         cmap='gray',
                         vmin=0,vmax=1);
            #ax_inh.set_yticks([],[])
            ax_inh.set_xticks([i+.5 for i,x in enumerate(self.names_inhibitors)])
            ax_inh.set_xticklabels(self.names_inhibitors, rotation=rotation,
                    fontsize=fcp)
            pylab.axis([0,self.inhibitors.shape[1], 0, self.inhibitors.shape[0]])

        # the stimuli
        if len(self.names_stimuli) > 0:
            pylab.sca(ax_stim_top)
            pylab.text(0.5,0.5, "Stimuli", color="blue", horizontalalignment="center",
                verticalalignment="center", fontsize=int(fcti/1.5))

        # the inhibitors
        if len(self.names_inhibitors) > 0:
            pylab.sca(ax_inh_top)
            pylab.text(0.5,0.5, "Inhibitors", color="blue",
                       horizontalalignment="center", verticalalignment="center",
                       fontsize=int(fcti/1.5))

        #colorbar
        # we build our own colorbar to place it on the RHS
        if colorbar and mode == "mse":
            pylab.sca(ax_cb)
            N = self._params['plot.colorbar.N']
            cbar = pylab.linspace(0, 1, N)
            indices = [int(x) for x in cbar**self.cmap_scale*(N-1)]

            cbar = [cbar[i] for i in indices]

            pylab.pcolor(np.array([cbar, cbar]).transpose(), cmap=cmap,
                         vmin=0, vmax=1);
            ax_cb.yaxis.tick_right()
            ticks = np.array(ax_cb.get_yticks())
            M = max(ticks)
            indices = [int(N*x) for x in ticks**self.cmap_scale/(M**self.cmap_scale)]
            ax_cb.set_yticks(indices)
            if vmax == 1:
                # set number of digits
                tic = np.array(indices)/float(N)
                tic = [int(x*100)/100. for x in tic]
                ax_cb.set_yticklabels(tic)
            else:
                ax_cb.set_yticklabels([int(x*100)/100. for x in np.array(indices)/float(N)*vmax])
            ax_cb.set_xticks([],[])

        pylab.sca(ax_main)
        return ax_main

    def _get_nExps(self):
         return len(self.experiments)
    nExps = property(_get_nExps, doc="return number of experiments")

    def _get_stimuli(self):
        return self.experiments.Stimuli
    stimuli = property(_get_stimuli, doc="return the stimuli dataframe")

    def _get_inhibitors(self):
        return self.experiments.Inhibitors
    inhibitors = property(_get_inhibitors, doc="return the inhibitors dataframe")

    def _get_nSignals(self):
        return len(self.df.columns)
    nSignals = property(_get_nSignals, doc="return the number of signals")

    def _get_xtlabels(self, logx=False):
        """build the time labels vector

        The vector is [t0,tmid,t0,tmid,...t0,tmid,tend]

        """
        t0 = self.times[0]
        t2 = self.times[-1]

        xtlabels = [int(t0),int((t2-t0)/2)] * self.nSignals + [int(t2)]
        return xtlabels

    def xplot(self,  bbox=False, *args, **kargs):
        """Same as :meth:`plot` using the xkcd layout for fun!"""
        with pylab.xkcd():
            self.plot(*args, **kargs)
            mybbox = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            if bbox:
                pylab.text(3, 12, "XMIDAS", fontsize=14,
                           verticalalignment='top', bbox=mybbox)

            tx = pylab.title("XMIDAS", bbox=dict(boxstyle='round',
                                                 facecolor='wheat', alpha=0.5))
            tx.set_position((0.6, 1.15))

    def plot(self, mode='data', **kargs):
        """Plot data contained in :attr:`experiment` and :attr:`df` dataframes.

        :param string mode: must be either "mse" or "data" (defaults to data)

        if mode is 'mse', calls also plot_layout, plot_data and plot_sim_data
        else calls plot_layout and plot_data

        .. plot::
            :include-source:
            :width: 80%

            from cno import XMIDAS, cnodata
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.plot(mode="mse")

        """
        pylab.clf()
        if mode == "mse":
            if self.df.min().min()<0:
                self.logging.warning("values are expected to be positive")
            if self.df.max().max()>1:
                self.logging.warning("values are expected to be normalised")
            kargs['mode'] = 'mse'
            self.plot_layout(**kargs)
            self.plot_data(**kargs)
            self.plot_sim_data(**kargs)
        else:
            kargs['mode'] = 'data'
            self.plot_layout(**kargs)
            self.plot_data(**kargs)

    def boxplot(self, mode="time"):
        """Plot boxplot of the dataframe with the data

        :param str mode: time or species

        .. plot::

            from cno import XMIDAS, cnodata
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.boxplot(mode="species")
            m.boxplot(mode="time")

        """
        if mode == "time":
             self.df.reset_index().boxplot(by="time")
             pylab.tight_layout()
        else:
            self.df.boxplot()

    def radviz(self, species=None, fontsize=10):
        """

        .. plot::
            :include-source:
            :width: 80%

            from cno import XMIDAS, cnodata
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.radviz(["ap1", "gsk3", "p38"])

        """
        if species == None:
            species = list(self.df.columns)
        from pandas.tools.plotting import radviz
        df = self.df.reset_index()
        del df['time']
        del df['cell']
        pylab.figure(1)
        pylab.clf()
        radviz(df[['experiment']+species], "experiment")
        pylab.legend(fontsize=fontsize)

    def to_midas(self, filename, expand_time_column=False):
        """Save XMIDAS into a MIDAS CSV file.

        :param str filename:

        """
        f = open(filename, "w")
        # TODO: cellline

        header = ["TR:%s:CellLine"%self.cellLine]

        try:
            header += ['TR:'+this for this in self.experiments.Stimuli.columns]
        except:
            pass
        try:
            header += ['TR:'+this+'i' for this in self.experiments.Inhibitors.columns]
        except:
            pass

        if expand_time_column == False:
            header += ["DA:ALL"]
        else:
            header += ["DA:" +  x for x in self.df.columns]

        try:
            header += ["DV:" +  x for x in self.df.columns]
        except Exception:
            print(self.df.columns)
            raise Exception

        header = ",".join(header)
        f.write(header + "\n")

        for time in self.times:
            # experiments of levels is 1
            # better to use experiments df so that order is same as in experiments
            for exp in self.experiments.index:

                #FIXME: if we drop an experiment, this fails. do we want to
                # update the laels when calling remove_experiment method
                measurements = self.df.xs((self.cellLine, exp, time))
                measurements = measurements[self.df.columns] # maybe not needed
                # keep it for now to be sure that order of measurements is same as in the header
                experiment = self.experiments.ix[exp]

                #FIXME use isinstance
                if type(measurements) == pd.core.series.Series:
                    if expand_time_column == False:
                        rowdata = [1] + list(experiment.values) + [time] + list(measurements)
                    else:
                        rowdata = [1] + list(experiment.values) +\
                            [time]*len(self.df.columns) + list(measurements)
                    f.write(",".join(["{}".format(x) for x in rowdata]) + "\n")

                elif type(measurements) == pd.core.frame.DataFrame:
                    for measurement in measurements.values:
                        if expand_time_column == False:
                            rowdata = [1] + list(experiment.values) + [time] + list(measurement)
                        else:
                            rowdata = [1] + list(experiment.values) +\
                                [time]*len(self.df.columns) + list(measurement)
                        f.write(",".join(["{}".format(x) for x in rowdata]) + "\n")
                else:
                    raise TypeError()

        f.close()

    def normalise(self, mode, inplace=True, changeThreshold=0, **kargs):
        """Normalise the data

        :param mode: 'time' or  'control'
        :param bool inplace:  Defaults to True.

        .. warning:: not fully tested. the mode "time" should work. The
            control mode has been tested on 2 MIDAS files only. This is
            a complex normalisation described in
            :class:`~cno.io.midas_normalisation.XMIDASNormalise`
        """
        if self.df.max().max() <=1  and self.df.min.min()>=0:
            raise CNOError("dataalready between 0 and 1")
        assert mode in ["time", "control"], "mode must be control or time"
        kargs['changeThreshold'] = changeThreshold
        if mode == "time":
            n = normalisation.XMIDASNormalise(self, **kargs)
            normed_midas = n.time_normalisation()
        elif mode == "control":
            n = normalisation.XMIDASNormalise(self, **kargs)
            normed_midas = n.control_normalisation()

        if inplace == False:
            return normed_midas
        else:
            self.df = normed_midas.copy()

    def to_obs(self):
        from cno.io.obs import OBS
        obs = OBS()

    def to_measurements(self):
        """Returns a Measurements instance

        Each datum in the dataframe :attr:`df` is converted into an
        instance of :class:`~cno.io.measurement.Measurements`.

        :return: list of experiments.

        ::

            mb = MIDASbuilder(m.to_measurements)
            mb.xmidas

        """
        experiments = []
        for row in self.df.iterrows():
            cellLine, exp, time = row[0]
            data = row[1]
            inhibitors = self.experiments.ix[exp].Inhibitors
            stimuli = self.experiments.ix[exp].Stimuli
            for species in data.index:
                e = Measurement(species, time=time, stimuli=dict(stimuli),
                        inhibitors=dict(inhibitors), value=data[species],
                        cellLine=cellLine)
                experiments.append(e)
        return experiments

    def add_uniform_distributed_noise(self, inplace=False, dynamic_range=1,
            mode="bounded"):
        """add random (uniformaly distributed) noise to the dataframe

        The noise is uniformly distributed between -0.5 and 0.5 and added to the
        values contained in the dataframe (for each combinaison of species and
        time/experiment). New values are :math:`\hat{X}=X + \mathcal{U}(-.5, .5)*dr`,
        where dr is the dynamical range. Note that final values may be below
        zero or above 1. If you do not want this feature, set the mode to
        "bounded" (The default is **free**). **bounded** means

        :param bool inplace: False by default
        :param float dynamic_range: a multiplicative value set to the noise
        :param str mode: bounded (between min and max of the current data) or free

        """
        # axis=1 means each values is modified
        # axis=0 means the same value is added to the entire column
        dr = dynamic_range
        # add a unique random number to each value irrespective of level/axis
        if mode == "bounded":
            # because we use min() and max(), we cannot use apply but must use
            # applymap. Otherwise, identical random are generated for each
            # species
            newdf = self.df.applymap(lambda x: x +
                    np.random.uniform(-x.min(),1-x.max())*dr)
        elif mode == "free":
            newdf = self.df + np.random.uniform(self.df) * dr
        else:
            raise CNOError("mode can be bounded or free")
        if inplace:
            self.df = newdf.copy()
        else:
            return newdf

    def add_gaussian_noise(self, mu=0, sigma=0.1, inplace=False):
        r"""add gaussian noise to the data. Results may be negative or above 1

        :param float beta: see equation
        :param float sigma: see equation (default to 0.1)
        :param bool inplace: Default to False

        .. math::

            \hat{X} = X + \mathcal{N}(\mu, \sigma)

        """
        # random.normal accepts x and takes its shape to return random values.
        # so, the addition of x and the ouptut of np.random.normal is as
        # expected: points by point.
        newdf = self.df.apply(lambda x:x + np.random.normal(mu, scale=sigma))
        if inplace:
            self.df = newdf
        else:
            return newdf

    def _make_df_compatible(self):
        """if you use MakeBuilder, you can really add any combi of experiments
        However, the resulting df built is not MIDAS compatible for sure. For example,
        not all same time are available for each experiment.

        """
        raise NotImplementedError

    def get_residual_errors(self, level="time", normed=False):
        r"""Return vector with residuals errors

        The residual errors are interesting to look at in the context
        of a boolean analysis. Indeed, residual errors is the minimum error
        that is unavoidable with a boolean network and comes from the discrete
        nature of such a model. In a boolean analysis, one would compare 0/1
        values to continuous values between 0 and 1. Therefore, however good
        is the optimisation, the value of the goodness of fit term cannot go
        under this residual error.

        :param: level to sum over
        :return: a time series with summed residual errors
            :math:`\sum (round(x)-x)^2` the summation is performed over species
            and experiment by default.

        .. doctest::

            >>> from cno import cnodata, XMIDAS
            >>> m = XMIDAS(cnodata("MD-ToyMMB_T2.csv"))
            >>> m.get_residual_errors()
            time
            0       0.000000
            10      2.768152
            100     0.954000
            dtype: float64

        if normed is False, returns:

        .. math::

            \sum_{level} \left( X - \bf{round}( X) \right)^2

        where level can be either over time or experiment.
        If normed is True, divides the time series by number of experiments
        and number of times

        .. note:: if normed set to False, same results as in
            CellNOptR for mode set to time.
        """
        diff = (self.df - self.df.apply(lambda x: x.round(), axis=1))
        diff_square = diff.apply(lambda x: x**2, axis=1)
        S = diff_square.sum(level=level).sum(axis=1)

        if normed :
            S /= len(self.experiments) * len(self.times)
        return S

    def copy(self):
        x = XMIDAS()
        x._data = self._data.copy()
        x._missing_time_zero = self._missing_time_zero
        x.df = self.df.copy()
        x._experiments = self.experiments.copy()
        x._cellLine = self.cellLine
        x.sim = self.sim.copy()
        x.errors = self.errors.copy()
        return x

    def __str__(self):
        txt = "Your data set contains {} cellLines\n".format(len(self.cellLines))
        txt += "   Current selected cell line is {}\n".format(self.cellLine)
        txt += "\nThe data contains \n"
        txt += "{} Species:\n\t {}\n".format(len(self.names_species),
                self.names_species)
        txt += "{} inhibitors:\n\t {}\n".format(len(self.names_inhibitors),
                self.names_inhibitors)
        txt += "{} stimuli:\n\t {}\n".format(len(self.names_stimuli),
                self.names_stimuli)
        return txt

    def pca(self, pca_components=2, fontsize=16):
        """PCA analysis

        .. plot::
            :include-source:
            :width: 80%

            from cno import XMIDAS, cnodata
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.pca()



        """
        from sklearn.decomposition import PCA
        pca = PCA(n_components=pca_components)
        #t = self.df[signal]
        #pca.fit(t.unstack(1))

        for signal in self.df.columns:
            t = self.df[signal].fillna(0)
            X = t.unstack()
            X_r = pca.fit(X).transform(X)
            pylab.plot(X_r[:,0], X_r[:,1], 'o', label=signal)
            print(signal, pca.explained_variance_ratio_)
        pylab.legend(fontsize=fontsize)

        return pca


    def discretize(self, **kargs):
        return self.discretise(**kargs)

    def discretise(self, inplace=True, N=2):
        """Discretise data by rounding up the values

        :param int N: number of discrete values (defaults to 2). If set to
            2, values will be either 0 or 1. If set to 5, values wil lbe in
            [0,0.25,0.5,0.75,1]
        :param inplace:

        .. warning. data has to be normalised
        """
        assert N >= 2, "N must be at least equal to 2"
        N = N-1.
        self.logging.info("Discretization between 0-1 assuming normalised data")
        if inplace:
            self.df = self.df.apply(lambda x: (x*N).round()/N)
        else:
            df = self.df.apply(lambda x: (x*N).round()/N)
            return df

    def round(self, inplace=True, decimals=0):
        """Round values to a given decimal"""
        if inplace == True:
            self.df = self.df.apply(lambda x : x.round(decimals=decimals))
        else:
            return self.df.values.apply(lambda x : x.round(decimals=decimals))

    def hcluster(self, mode="experiment", metric='euclidean', leaf_rotation=90,
                 leaf_font_size=12, **kargs):
        """Plot the hiearchical cluster (simple approach)

        .. plot::
            :include-source:
            :width: 80%

            from cno import XMIDAS, cnodata
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.hcluster("species")

        """
        kargs['leaf_font_size'] = leaf_font_size
        kargs['leaf_rotation'] = leaf_rotation

        assert mode in ["experiment", "time", "species"]
        from scipy.spatial.distance import pdist, squareform
        from scipy.cluster.hierarchy import linkage, dendrogram
        pylab.clf()
        if mode == "experiment":
            distxy = squareform(pdist(self.df.unstack("time"), metric=metric))
        elif mode == "time":
            distxy = squareform(pdist(self.df.unstack("experiment"), metric=metric))
        elif mode == "species":
            distxy = squareform(pdist(self.df.transpose(), metric=metric))

        try:
            R = dendrogram(linkage(distxy), **kargs)
        except:
            R = dendrogram(linkage(distxy, method='complete'), **kargs)

        self._debug = R
        leaves = R['leaves']
        if mode == "time":
            pylab.xticks(pylab.xticks()[0], self.times[leaves])
            pylab.title("Clustering by time")
        elif mode == "experiment":
            pylab.xticks(pylab.xticks()[0], self.experiments.index[leaves])
            pylab.title("Clustering by experiments")
        else:
            pylab.xticks(pylab.xticks()[0], self.species[leaves])
            pylab.title("Clustering by species")
        ylim = pylab.ylim()
        if ylim[0] == 0:
            pylab.ylim([0-ylim[1]*.05, ylim[1]])
        pylab.tight_layout()

    def heatmap(self, cmap="heat", transpose=False):
        """Hierarchical clustering on species and one of experiment/time level

        .. plot::
            :include-source:
            :width: 80%

            from cno import XMIDAS, cnodata
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.heatmap()

        .. note:: time zero is ignored. Indeed, data at time zero is mostly set
            to zero, which biases clustering.
        """
        #FIXME 1 looks like dendograms are not shown. why?
        from biokit.viz.heatmap import Heatmap
        data = self.df.query('time>0').unstack(2).ix[self.cellLine]
        if transpose:
            h = Heatmap(data.transpose())
        else:
            h = Heatmap(data)
        h.plot(cmap=cmap)
        return h

    def shuffle(self, mode="experiment", inplace=True):
        """Shuffle data

        this method does not alter the data but shuffle it around
        depending on the user choice.

        :param str mode: type of shuffling (see below)
        :param bool inplace: Defaults to True

        The **mode** parameter can be

        #. `timeseries` shuffles experiments and species; timeseries
            are unchanged.
        #. `all`: through times, experiments and species. No structure kept
        #. `signal` (or `species` or `column`): within a column, timeseries
           are shuffled. So, sum over signals is constant.
        #. `experiment` (or `index`): with a row (experiment),
           timeseries are shuffled.

        Original data:

        .. plot::
            :width: 80%

            from cno import XMIDAS, cnodata
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.plot()

        Shuffling all timeseries (shuffling rows and columsn in the plot):

        .. plot::
            :include-source:
            :width: 80%

            from cno import XMIDAS, cnodata
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.shuffle(mode="timeseries")
            m.plot()

        .. warning:: shuffling is made inplace.

        """
        if mode == "all":
            # The random.shuffle function does not work!! somehow sum of
            # data increases or decreases
            # One must use numpy.random.shuffle instead
            shape = self.df.shape
            data = self.df.values.reshape(shape[0]*shape[1])
            np.random.shuffle(data)
            count = 0
            # not very efficient but works for now
            for i in range(0, shape[0]):
                for j in range(0, shape[1]):
                    self.df.values[i][j] = data[count]
                    count += 1
        elif mode in ["signal", "species",  "column"]:
            # m.df.sum() is constant
            for c in self.df.columns:
                self.df[c] = np.random.permutation(self.df[c].values)
        elif mode in ["index", "experiment"]:
            # m.df.sum(level="experiment").sum(axis=1) is constant
            for this_index in self.df.index:
                np.random.shuffle(self.df.ix[this_index].values)
        elif mode == "timeseries":
            # swap boxes
            species = list(self.species)
            exps = list(self.experiments.index)
            pairs = []
            for s in species:
                for e in exps:
                    pairs.append((s,e))
            # permutation now
            np.random.shuffle(pairs)
            df = self.df.copy()
            # FIXME: cell not implemented yet
            count = 0
            for s in species:
                for e in exps:
                    rand_s = pairs[count][0]
                    rand_e = pairs[count][1]
                    data = list(df[rand_s].ix[self.cellLine].ix[rand_e])

                    self.df[s].ix[self.cellLine].ix[e] = data
                    count += 1
        else:
            raise CNOError("unknown mode {0} in shuffle()".format(mode))

    def sort_experiments_by_stimuli(self):
        """Sort experiments by stimuli

        Affects the experiment dataframe for th rendering but do not
        change the dataframe that contains the data.

        """
        N = len(self.experiments)
        list_exps = [('Stimuli', stimulus) for stimulus in self.names_stimuli]
        list_exps += [('Inhibitors', inhibitor) for inhibitor in self.names_inhibitors]
        new_order = self.experiments.groupby(list_exps).groups
        new_order = [new_order[x] for x in sorted(new_order.keys())]
        new_order = list(pylab.flatten(new_order))
        assert len(new_order )== N
        self._experiments = self.experiments.reindex_axis(new_order, axis=0)

    def sort_experiments_by_inhibitors(self):
        """Sort the experiments by inhibitors

        Affects the experiment dataframe for th rendering but do not
        change the dataframe that contains the data.
        """
        N = len(self.experiments)
        list_exps = [('Inhibitors', inhibitor) for inhibitor in self.names_inhibitors]
        list_exps += [('Stimuli', stimulus) for stimulus in self.names_stimuli]
        new_order = self.experiments.groupby(list_exps).groups
        new_order = [new_order[x] for x in sorted(new_order.keys())]
        new_order = list(pylab.flatten(new_order))
        assert len(new_order )== N
        self._experiments = self.experiments.reindex_axis(new_order, axis=0)

    def __eq__(self, other):
        if all(other.df == self.df) == False:
            return False
        if all(other.experiments == self.experiments) == False:
            return False
        return True


class Trend(object):
    """Utility that figures out the trend of a time series

    .. plot::
        :include-source:
        :width: 80%

        from cno import cnodata
        from cno.io import midas

        # get a time series
        xm = midas.XMIDAS(cnodata("MD-ToyPB.csv"))
        ts = xm.df['ap1']['Cell']['experiment_0']

        # set a trend instance
        trend = midas.Trend()
        trend.set(ts)
        trend.plot()

        trend.get_bestfit_color()


    """
    def __init__(self):
        self.data = None

    def set(self, ts):
        """Set the data with the parameter

        :param ts: a Pandas TimeSeries
        """
        self.data = ts

    def _get_times(self):
        return self.data.index.values
    times = property(_get_times, doc="return time array")

    def _get_values(self):
        return self.data.values
    values = property(_get_values, doc="return the value array")

    def _get_normed_times(self):
        # Here, time 0 is zero so we can divide by Max only
        return self.times / float(self.times.max())
    normed_times = property(_get_normed_times, doc="return normed time array")

    def _get_normed_values(self):
        # FIXME do we want to use the span ?? 
        #span = float(self.values.max()) - float(self.values.min())
        return self.values / float(self.values.max())
        #return (self.values -self.values.min())/ span
    normed_values = property(_get_normed_values, doc="return normed value array")

    def _get_alpha(self):
        try:
            alpha = np.trapz(self.normed_values, self.normed_times)
        except:
            alpha = 'white'
        return alpha
    alpha = property(_get_alpha, doc="return strength of the signal")

    def transient(self, x=None):
        if x == None:
            x = self.times
        M = max(x)
        v = np.array([(M-y)/(M/2.) if y>=M/2. else y/(M/2.) for y in x])
        return self._normed(v)

    def constant(self, x=0):
        v = np.array([x] * len(self.times))
        v = self._normed(v)
        return v

    def _normed(self, v):
        sumsq = np.sqrt(sum([this**2 for this in v]))
        return v/sumsq

    def earlier(self, x=None, n=3., N=4.):
        if x == None:
            x = self.times
        M = max(x)
        v = np.array([(M-y)/(n*M/N) if y>=M/N else y/(M/N) for y in x])
        return self._normed(v)

    def sustained(self, x=None, L=0.5):
        if x == None:
            x = self.times
        M = max(x)
        m = L * M
        v = np.array([y if y<m else m for y in x])
        return self._normed(v)

    def inverse_sustained(self, x=None, L=0.5):
        if x == None:
            x = self.times
        M = max(x)
        m = L * M
        v = np.array([(M-y) if y < m else M-m for y in x])
        return self._normed(v)

    def later(self, x=None, L=0.5):
        if x == None:
            x = self.times
        M = max(x)
        m = L * M
        v = np.array([0 if y<m else y-m for y in x])
        return self._normed(v)

    def _correlate(self, a, b):
        a = self._normed(a)
        b = self._normed(b)
        return np.correlate(a,b)

    def _get_correlation(self, a):
        correlation = {}
        correlation['later'] = self._correlate(a, self.later())
        correlation['earlier'] = self._correlate(a, self.earlier())
        correlation['earlier2'] = self._correlate(a, self.earlier(n=1,N=10))
        correlation['transient'] = self._correlate(a, self.transient())
        correlation['constant_half'] = self._correlate(a, self.constant(0.5))
        correlation['constant_unity'] = self._correlate(a, self.constant(1))
        correlation['sustained'] = self._correlate(a, self.sustained(L=.5))
        correlation['inverse_sustained'] = self._correlate(a, self.inverse_sustained(L=.5))
        return correlation

    def plot(self):
        """Plots the data and possible choices of trends"""
        data = self.values
        corrs = self._get_correlation(data)
        pylab.clf()
        pylab.plot(self.times, self._normed(data),
                   label="data", lw=4, ls="--", color='k')
        # transient
        pylab.plot(self.times, self.transient(), 'o-',
                   label="transient " + str(corrs['transient']))
        # earlier
        pylab.plot(self.times, self.earlier(), 'o-',
                   label="earlier " + str(corrs['earlier']))
        pylab.plot(self.times, self.earlier(n=1, N=10), 'o-',
                   label="earlier2 " + str(corrs['earlier2']))
        # later
        pylab.plot(self.times, self.later(), 'o-',
                   label="later " + str(corrs['later']))
        # constant
        pylab.plot(self.times, self.constant(.5), 'o-',
                   label="constant " + str(corrs['constant_half']))
        # sustained
        pylab.plot(self.times, self.sustained(L=.5), 'o-',
                   label="sustained" + str(corrs['sustained']))
        pylab.plot(self.times, self.inverse_sustained(L=.5), 'o-',
                   label="inv sustained" + str(corrs['inverse_sustained']))
        pylab.legend(fontsize=10)

    def get_bestfit(self):
        corrs = self._get_correlation(self.normed_values)
        keys,values = (corrs.keys(), corrs.values())
        #M  = max(values)
        return keys[np.argmax(values)]

    def get_bestfit_color(self):
        corrs = self._get_correlation(self.normed_values)
        keys,values = (corrs.keys(), corrs.values())
        #M  = max(values)
        res = keys[np.argmax(values)]

        if "constant" in res:
            return "black"
        elif "later" in res:
            return "red"
        elif "transient" in res:
            return "yellow"
        elif "earlier" in res:
            return "purple"
        elif "sustained" in res:
            return "green"
        else:
            return "white"


