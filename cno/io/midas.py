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

from easydev.logging_tools import Logging
import colormap

from midas_extra import Measurement
import midas_normalisation as normalisation

from cno.core import DevTools

__all__ = ["XMIDAS", "TypicalTimeSeries", 'MIDASReader']



class MIDAS(object):
    """This is s a description of the MIDAS format (let us call it 1.0)"""
    valid_codes = {
                   'ID':'identifier',
                   'TR':'stimuli/inhibitors',
                   'DA':'times',
                   'DV':'measurements'}

    ignore_codes = ['NOINHIB', 'NOCYTO', 'NOLIG', 'NO-CYTO', 'NO-INHIB', 'NO-LIG']

    def __init__(self, filename=None, verbose=False):
        """.. rubric:: Constructor

        :param str filename: the filename (including correct path)

        """
        self._header = []
        self._names = []
        self.verbose = verbose

        #: filename of the original data
        self.filename = filename

        self.cmap_scale = 0.5
        self.fontsize = 16
        self.logging = Logging("INFO")
        self._colormap = colormap.Colormap()

        self._dev = DevTools()


class MIDASReader(MIDAS):
    def __init__(self, filename, verbose='ERROR'):
        super(MIDASReader, self).__init__(filename, verbose)

    def read(self):
        if self.filename == None:
            return
        self._data = pd.read_csv(self.filename, skipinitialspace=True, sep=",")

        # figure out the cell line names
        self._preprocess_cellines()

        # remove columns that are invalid
        self._midas_validity()

        # some cleanup to remove columns that have to be ignored
        labels = ["TR:"+x for x in self.ignore_codes if "TR:"+x in self._data.columns]
        self._data = self._data.drop(labels, axis=1)

        #try:
        self._init()  # read MIDAS and convert to dataframe
        #except Exception, e:
        #    self.logging.warning("Could not interpret the MIDAS input data file")
        #    self.logging.warning(e.message)

    def _preprocess_cellines(self):
        #CellLine are tricky to handle with the MIDAS format because they use the
        #same prefix TR: as the treatments. You must be sure that (1) there are
        # 2 : signs (2) the suffix is called CellLine and (3) the middle name is
        # provided

        # Those could be recognised but only CellLine is correct. So, let us
        # enfore the correct one instead of dealing with all kind of user
        # choices.

        cellLines = [col for col in self._data.columns if "CellLine" in col]
        if len(cellLines) == 0:
            raise ValueError("Could not find any column with the required keyword 'CellLine'")

        for name in cellLines:
            if name.count(":") != 2:
                txt = "column name related to CellLine %s must have 2 : characters"
                txt += "{} has less or more than 2".format(name)
                raise ValueError(txt)

        # rename the cellLines if there are undefined
        # let us give thme the same name ('undefined')
        # this should be done at the level of the data
        for i,name in enumerate(cellLines):
            if name.split(":")[1] == "":
                self.logging.warning("Found a column related to CellLine without a name. Renamed to undefined.")
                columns = list(self._data.columns)
                columns[i] = "TR:undefined:CellLine"
                self._data.columns = columns

        # if there is only one undefined, no need for an error.
        # otherwise they will be an ambiguity.
        if len(cellLines) != len(set(cellLines)) and len(cellLines)>1:
            raise ValueError("some cellLines have the same name.")

        for this in self._data.columns:
            if this.split(":")[0].startswith("TR") == True:
                if "cellLine" in this.lower() and "CellLine" not in this:
                    raise ValueError("Found column with invalid tag. A Cell Line must be written 'TR:<name>:CellLine' not {}. Capitalisation matters.".format(this))

        celltype_names = self._get_cellLines()

        if len(cellLines) >=2:
            if self.cellLine == None or self.cellLine not in celltype_names:
                txt = "Error:: More than 1 celline was found.\n"
                txt += "You must select one amongst: {}".format([this.split(":")[1] for this in cellLines])
                raise ValueError(txt)
            else:
                # we could remove columns and rows where cell type is not correct.
                # but we are going to do it in the _init so that a user can
                # change his mind and reset the cell line to be looked at.
                pass

        if len(celltype_names) == 1:
            self._cellLine = celltype_names[0]

        #self._valid_cellLine_names = [this for this in self._]

    def _midas_validity(self):
        # checks are made of self._data only not df that will be built later on.

        if self._ignore_invalid_columns:
            for this in self._data.columns:
                columns = [this for this in self._data.columns if this[0:2] in self.valid_codes]
                bad = [this for this in self._data.columns if this[0:2] not in self.valid_codes]
                if len(bad):
                    self.logging.warning("Found columns that are invalid (do not start with {}). There are removed.".format(self.valid_codes.keys()))
                    self._data = self._data[columns]

            #columns = [this in self._data.columns if this.startswith("")]
            pass

        # check validity of TR:
        for this in self._data.columns:
            if ":" not in this:
                txt = "Error in header of the input MIDAS file\n"
                txt += "    Column's name must contain the special character ':'\n"
                raise ValueError(txt)
            if this.split(":")[0] not in self.valid_codes.keys():
                txt = "Error in header of the input MIDAS file\n"
                txt += "    Column's name must start with one of the valid code {}\n".format(self.valid_codes)
                raise ValueError(txt)


        # check if zero time data is available otherwise need to call
        #maybe better to play with -df
        #_duplicate_time_zero_using_inhibitors_only
        df = self._data[[this for this in self._data.columns if this.startswith("DA")]]
        unique_times = list(set(df.as_matrix().flatten()))
        unique_times.sort()
        if len(unique_times) <2:
            raise ValueError("Must contains at least 2 time points including time zero")
        if 0 not in unique_times or len(unique_times)<=1:
            raise ValueError("You must have zero times in the MIDAS file, that was not found.")
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
            raise ValueError("Header of MIDAS file has no columns starting with DV. expects at least one")


class XMIDAS(MIDASReader):
    """XMIDAS dat structure. X stands for extended and replaces
    :class:`MIDASReader` class.

    ::

        from cellnopt.core import XMIDAS
        m = XMIDAS(cnodata("tes.csv"))
        m.df # access to the data frame
        m.scale_max(gain=0.9)   # scale over all experiments and time finding the max
                        # and scaling (divide by max) for each species individually.
        m.corr()        # alias to m.df.corr() removing times/experiments
        columns



        tuples = [(exp, time) for exp in ["exp1", "exp2"] for time in [0,1,2,3,4]]
        index = pd.MultiIndex.from_tuples(tuples, names=["experiment", "time"])
        xx = pd.DataFrame(randn(10,2), index=index, columns=["Akt", "Erk"])


    What remains to be done ?

        * average over celltypes  To be done in MultiMIDAS


    .. warning:: if there are replicates, call average_replicates before creating a
        simulation or calling plot"mode="mse")

    .. todo:: when using MSE, an option could be to average the errors by taking into
        account the time. In other words, a weight/integral. if t = 1,2,3,4,5,10,60, the
        errors on 1,2,3,4,5 are more  important than between 5,10,60. does it make sense ?

    .. todo:: make df a property to handle sim properly, if not scaled,
        sim and exp seems to have the same scale also errors are large as expected

    .. warning:: MD-TR-33333-JITcellData.csv contains extra ,,,, at the end. should
        be removed or ignored

    .. todo:: colorbar issue with  midas.XMIDAS("share/data/MD-test_4andgates.csv")


    .. todo:: when plotting, if there is only 1 stimuli and 5-6 inhibitors, the
        width of the stimuli is the same as the one with the inhibitors. cell size
        should be identical, not stretched. See e..g., "EGFR-ErbB_PCB2009.csv"
:

    .. todo:: when ploting the mse, we should be able to plotonly a subset of the time
      indices (useful for bollean analysis at a given time)


    .. todo::MIDAS have two ways of coding stimuli/inhibitors a short and long version
        TR:aa / TR::aai or TR:aa:Stimuli / TR:aa:Inhibitors note that in the shotr case,
        the letter i is used to encode inhibitor, which is not robust at all.


    ..todo:: a MIDAS class to check validity just to simplfy the XMIDAs class itself.

    .. todo:: inhibitors ends in :i to avoid clashes with same name in stimuli..


    API:

        - you can store several cell line within one MIDAS file. However, XMIDAS
          handles only 1 for now.

    """
    def __init__(self, filename=None, cellLine=None, verbose=False):
        super(XMIDAS, self).__init__(filename)

        self._cellLine = cellLine

        # multi index related
        # position of the columns for the multi index
        self._celltype_index = 0
        self._experiment_index = 1
        self._time_index = 2
        # names of the multi index level (rows)
        self._levels = ["cell", "experiment", "time"]

        self.verbose = verbose

        self._ignore_invalid_columns = True


        self._data = pd.DataFrame()
        self._experiments = pd.DataFrame()
        self.df = pd.DataFrame()

        self.read()

        self.create_empty_simulation()
        self.errors = self.sim.copy()
        self._missing_time_zero = False

    def _manage_replicates(self):
        """

         tuples = [("HepG2", exp, time) for exp in ["exp1", "exp2"] for time in   [0,1,2,3,4]]
         tuples += [("Liver", exp, time) for exp in ["exp1", "exp2"] for time in [0,1,2,3,4]]
         tuples += [("HepG2", exp, time) for exp in ["exp1"] for time  in [0,1,2,3,4]]
         index = pd.MultiIndex.from_tuples(tuples, names=["CellType", "experiment", "time"])
         xx = pd.DataFrame(randn(25,2), index=index, columns=["Akt", "Erk"])

        xxx = xx.groupby(level=["CellType", "experiment", "time"]).agg([mean, std])
        xxx.to_csv("test.csv", tupleize_cols=False)
        pd.read_csv("test.csv", index_col=[0,1,2], tupleize_cols=False, header=[0,1])

        x.df = x.df[[this for this in x.df.columns if "mean" in this]]
        x.df.columns = [this[0] for this in x.df.columns]

        """
        groups = self.df.groupby(level=self._levels).groups
        if any([len(this)>1 for this in groups.values()])==False:
            self.logging.info("No replicates found")
        else:
            newdf = self.df.groupby(level=self._levels).agg([np.mean, np.std])
            return newdf

    def average_replicates(self, inplace=False):
        df = self._manage_replicates()

        if isinstance(df,pd.DataFrame):
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
        """

        .. doctest::

            >>> m = XMIDAS("MD-ToyPB.csv")
            >>> m['p38']
            >>> m['Cell','experiment_0', '0']

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

    def _get_cues(self):
        #cues = [x for x in self.experiments.columns if x.startswith('TR:')]
        cues = [x for x in self.experiments.columns]
        return cues
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
    cellLines = property(_get_cellLines)

    def _get_cellLine(self):
        return self._cellLine
    def _set_cellLine(self, name):
        # TODO check valid name
        names = self.cellLines
        if name not in names:
            raise ValueError("Invalid cellLine name {}. Valid ones are {}".format(name, names))
        self._cellLine = name
        # TODO: do we need to call _init again ?
        self._init()
    cellLine = property(_get_cellLine, _set_cellLine)

    def _get_times(self):
        times = self.df.index.levels[self._time_index]
        return sorted(list(times))
    times = property(_get_times)

    def _get_names_inhibitors(self):
        try:
            return list(self.experiments.Inhibitors.columns)
        except:
            return []
    names_inhibitors = property(_get_names_inhibitors)

    def _get_names_stimuli(self):
        try:
            return list(self.experiments.Stimuli.columns)
        except:
            return []
    names_stimuli = property(_get_names_stimuli)

    def _init(self):

        # select only data that matches the cell line choice made by the user.
        cellLine = 'TR:%s:CellLine' % self.cellLine
        _data = self._data[self._data[cellLine] == 1]
        # and remove all column with the CellLine keyword
        _data = _data[ [col for col in _data.columns if "CellLine" not in col]]
        #drop ID columns if any
        _data = _data[ [col for col in _data.columns if col.startswith("ID")==False]]


        df_tr = _data[[this for this in _data.columns if this.startswith("TR")]]
        df_da = _data[[this for this in _data.columns if this.startswith("DA")]]
        df_dv = _data[[this for this in _data.columns if this.startswith("DV")]]
        # TODO sort alphabetical ignoring big caps
        df_dv = _data[[this for this in _data.columns if this.startswith("DV")]]

        value_experiments = _data[df_tr.columns].copy()
        value_experiments.replace("NaN", 0, inplace=True)

        value_signals = _data[df_dv.columns].as_matrix()
        value_times = _data[df_da.columns]

        names = [this for this in df_tr[:] if "CellLine" not in this]
        self._experiments = _data[names].drop_duplicates()
        self._experiments.index = range(0, self._experiments.shape[0])
        self._experiments.index = ["experiment_{}".format(this) for this in  self._experiments.index]

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
        self._experiments = self._experiments.applymap(lambda x: 0 if isinstance(x, basestring) and x.isspace() else x)
        self._experiments = self._experiments.convert_objects(convert_numeric=True, copy=True)

        index =  pd.MultiIndex.from_tuples(tuples, names=self._levels)
        #keep = [this for this in self.df.columns if this not in ["experiments", "times"]]
        names_species = [x for x in _data.columns if x.startswith('DV:')]
        names_species = [x[3:] for x in names_species]

        self.df = pd.DataFrame(value_signals, index=index, columns=names_species)
        self.df = self.df.sortlevel(["experiment"])
        self.df = self.df.sort_index(axis=1) # sort the species

        if self._missing_time_zero == True:
            self._duplicate_time_zero_using_inhibitors_only()

        if self.df.shape[0] > len(self.times) * self.experiments.shape[0]:
            self.logging.warning("WARNING:: you may have duplicated experiments, pleiase average the replaicates using self.average_replicates(inplace=True)")

        if self.df.max(skipna=True).max(skipna=True) > 1:
            self.logging.warning("WARNING:: values larger than 1. You may want to normalise/scale the data")

        # Get rid of TR in experiments
        self._experiments.columns = [this.replace("TR:", "") for this in self._experiments.columns]
        #

        cues = []
        for c in self._experiments.columns:
            if c.count(":") >=2:
                raise ValueError("Invalid header. Found more than 2 : sign in a column")

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
        self._experiments = pd.DataFrame(self._experiments.values,
                     columns=[['Stimuli']*len(stim) + ['Inhibitors']*len(inh),
                              [x.replace(":i","") for x in self._experiments.columns]],
                        index= self._experiments.index)

    def _check_consistency_data(self):
        # times consistency all times must have same length
        #
        pass

    def _duplicate_time_zero_using_inhibitors_only(self):
        """
        Sometimes the time zero data sets are not explicitly written in MIDAS
        files. One example is MD-ExtLiverHepG2-MCP2010-mod4.csv from the

        """
        self.logging.warning("WARNING:: duplicating time zeros data to fit experiment at other times")
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
                these_inhibitors =  self.experiments.ix[this_exp][self.inhibitors]

                for this_exp_intern in experiments:
                    times = self.df.ix[self.cellLine].ix[this_exp_intern].index
                    if 0 in times:
                        if all(self.experiments.ix[this_exp_intern][these_inhibitors.index] ==
                            these_inhibitors):
                            #print("{}(no time zero found) is similar to {}".format(this_exp,  this_exp_intern))
                            break # so that if there are several replicates found, we pick up the first one only

                # get the times for this experimenti. it must contain the time zero
                newdata = self.df.xs((self.cellLine, this_exp_intern))
                # we only need the time 0
                newrow = newdata[newdata.index == 0]
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
        df = df.sortlevel(self._levels[1]) # experiment
        self.df = df.copy()

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
            self.sim.set_index(self.levels, inplace=True)
            self.errors.set_index(self._levels, inplace=True)

    def remove_cellLine(self, labels):
        """Remove a cellLine from the dataframe.

        Does not really work since there is only one cellLine in the dataframe.
        all data is contained in :attr:`data` but the current dataframe contains only
        one, which can be changed simply by setting the cellLine attribute with one of the
        valid cellLine found in the :attr:`cellLines` attribute

        """
        self._remove_labels_from_level(labels, self._levels[0])

    def remove_times(self, labels):
        """Remove time values from the data

        :param list labels: one time point or a list of time points. Valid time points
            are in the :attr:`times` attribute.

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
        labels = self._str_or_list_to_list(labels)
        if level == "experiment":
            labels = [x if "experiment" in str(x) else "experiment_"+str(x)
                    for x in labels]

        self.reset_index()
        for label in labels:
            if label not in set(self.df[level]):
                self.logging.warning("{} not in times. Skipped".format(label))
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
        self._remove_column_experiment(labels)

    def remove_inhibitors(self, labels):
        """Remove inhibitor(s) from the :attr:`experiment` dataframe

        :param labels: a string or list of string representing the inhibitor(s)
        """
        self._remove_column_experiment(labels)

    def _remove_column_experiment(self, labels):
        labels = self._str_or_list_to_list(labels)
        for label in labels:
            if label not in set(self.experiments.columns):
                self.logging.warning("{} not in times. Skipped".format(label))
            else:
                self.experiments.drop(labels, axis=1, inplace=True)

    def rename_stimuli(self, names_dict):
        """Rename stimuli in the :attr:`experiment` dataframe

        :param names_dict: a dictionary with names (keys) to be replaced (values)

        ::

            from cellnopt.core import *
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.rename_species({"erk":"ERK", "akt":"AKT"})

        .. seealso:: :meth:`rename_stimuli`, :meth:`rename_inhibitors`

        """
        self._dev.check_param_in_list(names_dict.keys(), self.experiments.columns)
        columns = list(self.experiments.columns)
        columns = [c if c not in names_dict.keys() else names_dict[c] for c in columns]
        self.experiments.columns = columns

    def rename_inhibitors(self, names_dict):
        """Rename inhibitors

        :param names_dict: a dictionary with names (keys) to be replaced (values)

        ::

            from cellnopt.core import *
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.rename_species({"raf:i":"RAF:i"})

        .. seealso:: :meth:`rename_stimuli`, :meth:`rename_species`

        .. warning:: inhibitor name must end with the string **:i**

        .. todo:: sanity check that the pair of key/value contain the :i characters


        """
        self._dev.check_param_in_list(names_dict.keys(), self.experiments.columnss)
        columns = list(self.experiments.columns)
        columns = [c if c not in names_dict.keys() else names_dict[c] for c in columns]
        self.experiments.columns = columns

    def rename_species(self, names_dict):
        """Rename species in the main :attr:`df` dataframe

        :param names_dict: a dictionary with names (keys) to be replaced (values)

        ::

            from cellnopt.core import *
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.rename_species({"erk":"ERK", "akt":"AKT"})


        .. seealso:: :meth:`rename_stimuli`, :meth:`rename_inhibitors`


        """
        self.dev.check_param_in_list(names_dict.keys(), self.df.columns)
        columns = list(self.df.columns)
        columns = [c if c not in names_dict.keys() else names_dict[c] for c in columns]
        self.df.columns = columns
        self.sim.columns = columns

    def rename_cellLine(self, to_replace):
        """Rename cellLine indices

        :param dict to_replace: dictionary with mapping of values to be replaced.

        For example; to convert time in minutes to time in seconds, use something like::

            m.rename_cellLine({"undefined": "PriHu"})

        """
        self.reset_index()
        self.df.replace({self._levels[0]: to_replace}, inplace=True)
        self.set_index()

    def rename_time(self, to_replace):
        """Rename time indices

        :param dict to_replace: dictionary with mapping of values to be replaced.

        For example; to convert time in minutes to time in seconds, use something like::

            m.rename_time({0:0,1:1*60,5:5*60})

        """
        self.reset_index()
        self.df.replace({"time": to_replace}, inplace=True)
        self.set_index()

    def merge_times(self, how="mean"):
        raise NotImplementedError

    def add_experiment(self, e):
        raise NotImplementedError

    def corr(self, names=None, cmap=None):
        """plot correlation between the measured species

        :param list names: restriction to some species if provided.
        :param string cmap: a valid colormap (e.g. jet). Can also use "green" or "heat".

        .. plot::
            :include-source:
            :width: 80%

            >>> from cellnopt.core import *
            >>> m = XMIDAS(cnodata("MD-ToyPB.csv"))
            >>> m.corr(cmap="green")

        """
        cmap = self._get_cmap(cmap)

        corr = self.df.corr()
        N = corr.shape[0]
        names = self.df.columns[:]

        pylab.clf()
        pylab.pcolor(corr, edgecolors="k", cmap=cmap);
        pylab.xticks([x+0.5 for x in range(0,N)], names, rotation=90)
        pylab.yticks([x+0.5 for x in range(0,N)], names)
        pylab.tight_layout()
        pylab.colorbar()
        return corr

    def scale_max_across_experiments(self, inplace=True):
        """Divide each species column by max acrosss all experiments

        In the MIDAS plot, this is equivalent to dividig each column by
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

            X = \frac{X-m}{M-m}


        """
        m = self.df.min()
        M = self.df.max()
        data = (self.df - m)/(M-m.astype(np.float64))
        if inplace:
            self.df = data.copy()
        else:
            return data

    def scale_max_by_experiments(self, inplace=True):
         newdf = self.df.divide(self.df.max(level="experiment"), level="experiment")
         if inplace:
             self.df = newdf.copy()
         else:
             return newdf

    def scale_min_max_by_experiments(self, inplace=True):
        m = self.df.min(level="experiment")
        M = self.df.max(level="experiment")
        newdf = self.df.sub(m, level="experiment")
        newdf = newdf.divide(M-m, level="experiment")
        if inplace:
            self.df = newdf.copy()
        else:
            return newdf

    def scale_max(self, inplace=True):
        """Divide all data by the maximum over entire data set"""
        M = self.df.max().max()
        if inplace:
            self.df /= M
        else:
            return self.df / M

    def scale_min_max(self, inplace=True):
        r"""Divide all data by the maximum over entire data set

        .. math::

            X = \frac{X-m}{M-m}

        where :math:`m = min_{e,s,t} X` and :math:`M = max_{e,s,t} X`,
        with :math:`e` the experiment, with :math:`s` the species,
        with :math:`t` the time.

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
        self.sim = self.df *0 + numpy.random.uniform(size=self.df.shape)

    def get_diff(self, sim=None, norm="square", normed=True):
        """

        return difference between X and simulation. Take absolute.
        if norm == square or norm == absolute takes absolute values.
        if norm == square, also take power of 2.

        divide by number of time points.

        .. todo:: doc

        """
        # dataframe cannot be compared to None so, we need this trick:
        assert norm in ["absolute", "square"]
        if isinstance(sim, types.NoneType):
            sim = self.sim

        if norm == "square":
            diff = (sim - self.df).abs()**2
        else:
            diff = (sim - self.df).abs()

        diff = diff.sum(level="experiment")

        if normed:
            N = len(self.times)
            diff = diff/float(N)
        return diff

    def plotSim(self, markersize=3, logx=False, linestyle="--", lw=1,
            color="b", marker="x", **kargs):
        """plot experimental curves

        .. plot::
            :width: 80%
            :include-source:

            >>> from cellnopt.core import *
            >>> m = midas.MIDASReader(cnodata("MD-ToyPB.csv"));
            >>> m.plotMSEs()
            >>> m.plotExp()
            >>> m.plotSim()

        """
        times = numpy.array(self.times)
        # if simulation do not have the same number of points as data
        simtimes = numpy.array(self.sim.index.levels[2])

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

    def plotExp(self, markersize=3, logx=False,color="black", **kargs):
        """plot experimental curves

        .. plot::
            :width: 80%
            :include-source:

            >>> from cellnopt.core import *
            >>> m = midas.MIDASReader(cnodata("MD-ToyPB.csv"));
            >>> m.plotMSEs()
            >>> m.plotExp()

        .. note:: called by :meth:`plot`
        .. seealso:: :meth:`plot`, :meth:`plotMSEs`, :meth:`plotSim`
        """
        mode = kargs.get("mode", "trend")
        normalise = kargs.get("normalise", True)
        times = np.array(self.times)
        max_time = float(max(self.times))

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

        #for isim, sim in enumerate(self.sim):

        # vMax over all data
        if normalise:
            vMax = float(self.df.max(skipna=True).max(skipna=True))
        else:
            vMax = 1.

        #norm = np.trapz([1]*len(self.times), self.times/max_time)
        #print(norm)

        if mode == "trend":
            ts = TypicalTimeSeries(self.times)

        # TODO must be using the index instead of a range ince indices may not
        # start at zero
        for i in range(0, len(self.experiments)):
            #vMax = float(self.df.max(skipna=True).max(skipna=True))
            for j in range(0, self.nSignals):
                # divide data by 1.1 to see results close to 1.

                y = self.df[self.names_species[j]][self.cellLine][self.experiments.index[i]]
                times = numpy.array(y.index) / max_time
                if mode == "trend":
                    ts._times = times
                #vMax = self.df[self.names_species[j]][self.cellLine].max()

                y = y / vMax / 1.05

                #y = numpy.array([x[i,j] for x in self.exp])/vMax/1.05
                if mode == "trend":
                    try:
                        # time must be normalised by max so alpha is <=1
                        alpha = np.trapz(y.values, times)
                    except:
                        color="white"

                    try:
                        color = ts.get_bestfit_color(y.values)
                    except Exception:
                        color="white"

                    if color == "white":
                        colorc = "k"
                    else:
                        colorc=color

                    try:
                        pylab.plot(times+j, y+self.nExps-i-1 , 'k-o',
                                   markersize=markersize, color=colorc)
                        pylab.fill_between(times+j, y+self.nExps-1-i ,
                                           self.nExps-1-i, alpha=alpha,
                                           color=color)
                    except:
                        pass
                else:
                    pylab.plot(times+j, y+self.nExps-i-1 , 'k-o',
                               markersize=markersize, color="k")

                #    plot(times+j, sim[i,j]/1.05+(self.nExps-i-1), 'b--o', markersize=markersize)
        pylab.gca().set_xticklabels(xtlabels, fontsize=kargs.get("fontsize",10))
        pylab.gca().set_xticks(xt)

    def _get_cmap(self, cmap=None):
        if cmap == "heat":
            cmap = self._colormap.get_cmap_heat_r()
        elif cmap == "green":
            cmap = self._colormap.get_cmap_red_green()
        return cmap

    def plotMSEs(self, cmap="heat", N=10, norm="square",
        rotation=90,margin=0.05, colorbar=True, vmax=None, vmin=0.,
        mode="trend", **kargs):
        """plot MSE errors and layout


        .. plot::
            :width: 80%
            :include-source:

            >>> from cellnopt.core import *
            >>> m = midas.MIDASReader(cnodata("MD-ToyPB.csv"));
            >>> m.plotMSEs()

        .. todo:: error bars

        .. todo:: dynamic fontsize in the signal names ?

        .. note:: called by :meth:`plot`
        .. seealso:: :meth:`plot`, :meth:`plotMSEs`, :meth:`plotSim`

        .. todo:: need to make it more modular e.g. no cues matrices
        """
        if mode == "trend":
            #should be one with zero being white
            cmap = self._colormap.get_cmap_heat_r()
        else:
            cmap = self._get_cmap(cmap)


        diffs = self.get_diff(self.sim, norm=norm)
        diffs = diffs.ix[self.experiments.index]

        pylab.clf();

        bW = 0.1
        cH = 0.1

        if len(self.names_inhibitors)>0:
            bbW = 0.1
        else:
            bbW = 0

        aH = 1-cH-4*margin
        aW = 1-bW-5*margin - bbW

        # MAIN subplot with signals
        a = pylab.axes([margin, 2*margin, aW, aH])
        M = numpy.nanmax(diffs) # figure out the maximum individual MSE
        m = numpy.nanmin(diffs) # figure out the minimum individual MSE
        vmax_user = vmax
        vmax= max(1, M)         # if M below 1, set the max to 1 otherwise to M
        if vmax_user:
            vmax = vmax_user

        if mode == "mse":
            diffs = masked_array = np.ma.array (diffs, mask=np.isnan(diffs))
            cmap.set_bad("grey", 1.)
            pylab.pcolormesh(pylab.flipud(diffs)**self.cmap_scale, cmap=cmap, vmin=vmin, vmax=vmax, edgecolors='k');
        elif mode == "trend":
            cmap.set_bad("grey", 1.)
            pylab.pcolor(pylab.flipud(diffs*0), cmap=cmap, edgecolors='k');

        a.set_yticks([],[])
        pylab.axis([0, diffs.shape[1], 0, diffs.shape[0]])
        # Could add the names
        ax2 = a.twiny()
        ax2.set_xticks([i+.5 for i,x in enumerate(self.names_species)])
        N = len(self.names_species)
        ax2.set_xticks(pylab.linspace(0.5,N-1, N))
        ax2.set_xticklabels(self.names_species, rotation=90)

        # the stimuli
        b = pylab.axes([margin*2+aW, 2*margin, bW, aH])
        stimuli = numpy.where(numpy.isnan(self.stimuli)==False, self.stimuli, 0.5)

        pylab.pcolor(1-pylab.flipud(stimuli), edgecolors='gray', cmap='gray',vmin=0,vmax=1);
        b.set_yticks([],[])
        b.set_xticks([i+.5 for i,x in enumerate(self.names_stimuli)])
        b.set_xticklabels(self.names_stimuli, rotation=rotation)
        pylab.axis([0,self.stimuli.shape[1], 0, self.stimuli.shape[0]])

        # the inhibitors
        if len(self.names_inhibitors)>0:
            bb = pylab.axes([margin*5+aW, 2*margin, bbW, aH])
            inhibitors = numpy.where(numpy.isnan(self.inhibitors)==False, self.inhibitors, 0.5)
            pylab.pcolor(1-pylab.flipud(inhibitors), edgecolors='gray', cmap='gray',vmin=0,vmax=1);
            bb.set_yticks([],[])
            bb.set_xticks([i+.5 for i,x in enumerate(self.names_inhibitors)])
            bb.set_xticklabels(self.names_inhibitors, rotation=rotation)
            pylab.axis([0,self.inhibitors.shape[1], 0, self.inhibitors.shape[0]])


        d = pylab.axes([margin*2+aW, margin*3+aH, bW, cH])
        pylab.text(0.5,0.5, "Stimuli", color="blue", horizontalalignment="center",
            verticalalignment="center", fontsize=self.fontsize)
        #pcolor(1-numpy.zeros((1, 1)), edgecolors='b', cmap='gray', vmax=1, vmin=0);
        d.set_xticks([],[])
        d.set_yticks([],[])

        if len(self.names_inhibitors)>0:
            dd = pylab.axes([margin*5+aW, margin*3+aH, bbW, cH])
            pylab.text(0.5,0.5, "Inhibitors", color="blue", horizontalalignment="center",
                verticalalignment="center", fontsize=self.fontsize)
            #pcolor(1-numpy.zeros((1, 1)), edgecolors='b', cmap='gray', vmax=1, vmin=0);
            dd.set_xticks([],[])
            dd.set_yticks([],[])

        #colorbar
        # we build our own colorbar to place it on the RHS
        if colorbar and mode=="mse":
            e = pylab.axes([margin*3.5+aW+bW+bbW, 2*margin, margin/2, aH])
            cbar = pylab.linspace(0, 1, N)
            indices = [int(x) for x in cbar**self.cmap_scale*(N-1)]

            cbar = [cbar[i] for i in indices]

            pylab.pcolor(numpy.array([cbar, cbar]).transpose(), cmap=cmap, vmin=0, vmax=1);
            #d.set_xticks([],[])
            e.yaxis.tick_right()
            #e.yaxis.xticks([0,1][0,1])
            # todo: why is it normalised by 20?


            ticks = numpy.array(e.get_yticks())
            M = max(ticks)
            indices = [int(N*x) for x in ticks**self.cmap_scale/(M**self.cmap_scale)]
            e.set_yticks(indices)
            if vmax == 1:
                # set number of digits
                tic = numpy.array(indices)/float(N)
                tic = [int(x*100)/100. for x in tic]
                e.set_yticklabels(tic)
            else:
                e.set_yticklabels([int(x*100)/100. for x in numpy.array(indices)/float(N)*vmax])

            e.set_xticks([],[])

        pylab.sca(a)

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
        """Same as :meth:`plot` using the xkcd layout !"""
        with pylab.xkcd():
            self.plot(*args, **kargs)
            mybbox = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            if bbox:
                pylab.text(3, 12, "XMIDAS", fontsize=14,
                           verticalalignment='top', bbox=mybbox)

            tx = pylab.title("XMIDAS", bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            tx.set_position((0.6, 1.15))

    def plot(self, **kargs):
        """

        :param string mode: must be either "mse" or "trend" (defaults to trend)

        calls plotMSEs and plotExp

        if mode == mse, calls also plotSim

        .. plot::
            :include-source:
            :width: 80%

            from cellnopt.core import *
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.plot(mode="trend")

        .. todo:: a zero line

        .. no stimuli or inhibitors
        """
        mode = kargs.get("mode", "trend")
        kargs['mode'] = mode
        assert mode in ["mse", "trend"]

        if mode == "mse":
            if self.df.min().min()<0:
                self.logging.warning("values are expected to be positive")
            if self.df.max().max()>1:
                self.logging.warning("values are expected to be normalised")
        self.plotMSEs(**kargs)

        self.plotExp(**kargs)
        if mode == "mse":
            self.plotSim(**kargs)

    def save2midas(self, filename, expand_time_column=False):
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

    def _todelete_get_control(self, experiment_name):
        inhibitors = self.inhibitors.ix[experiment_name]
        # here is a sub selction where experiment matches the experiment_name
        mask = (self.inhibitors == inhibitors).all(axis=1)
        # indices contains all experiment that have the same inhibitors
        indices = self.inhibitors[mask].index

        # now from those experiment, which one is the control (i.e., all stimuli are off)
        stimuli = self.stimuli.ix[indices]
        mask = stimuli.sum(axis=1)==0
        control = stimuli[mask].index

        # TODO assert control is unique
        assert len(control) == 1
        return control[0]

    def normalise(self, mode, inplace=True, changeThreshold=0, **kargs):
        """Normalise the data

        :param mode: time or  controle
        :param bool inplace:  Defaults to True.

        see :mod:`normalise.XMIDASNormalise`

        .. warning:: not fully tested. the mode "time" should work. The
            control mode has been tested on 2 MIDAS file only.
        """
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

    def to_measurements(self):
        """Returns a Measurements instance
        
        Each datum in the dataframe :attr:`df` is converted into an
        instance of :class:`~cellnopt.core.xmidas.Experiment`.

        :return: list of experiments.

        ::

            mb = MIDASbuilder(m.to_measurements)
            mb.xmidas

        """
        experiments = []
        for row in self.df.iterrows():
            cellLine, exp, time = row[0]
            data = row[1]
            inhibitors = self.experiments.ix[exp][self.inhibitors]
            stimuli = self.experiments.ix[exp][self.stimuli]
            for species in data.index:
                e = Measurement(species, time=time, stimuli=dict(stimuli),
                        inhibitors=dict(inhibitors), measurement=data[species],
                        cellLine=cellLine)
                experiments.append(e)
        return experiments

    def add_uniform_distributed_noise(self, inplace=False, dynamic_range=1,
            mode="bounded"):
        """add random (uniformaly distributed) noise to the dataframe

        The noise is uniformliy distributed between -0.5 and 0.5 and added to the
        values contained in the dataframe (for each combinaison of species and
        time/experiment). New values are :math:`\hat{X}=X + noise(-.5, .5)*dr`,
        where dr is the dynamical range. Note that final values may be below
        zero or above 1. If you do not want this feature, set the mode to
        "bounded" (The default is **free**). **bounded** mens

        :param bool inplace: False by default
        :param float dynamic_range: a multiplicative value set to the noise
        :param bool min_value: final values below min are set to min (default is 0)
        :param bool max_value: final values above max are set to max (default is 1)

        """
        # axis=1 means each values is modified
        # axis=0 means the same value is added to the entire column
        dr = dynamic_range
        # add a unique random number to each value irrespective of level/axis
        # N = len(self.df.index)
        # x.df = x.df.apply(lambda x: x + np.random.normal(0, 10, size=N), axis=0)
        # OR change N to be column length and loop over axis 1
        # N = len(self.df.columns)
        # x.df = x.df.apply(lambda x: x + np.random.normal(0, 10, size=N), axis=1)
        if mode=="bounded":
            # because we use min() and max(), we cannot use apply but must use
            # applymap. Otherwise, identical random are generated for each
            # species
            newdf = self.df.applymap(lambda x: x +
                    np.random.uniform(-x.min(),1-x.max())*dr)
        elif mode == "free":
            newdf = self.df + np.random.uniform(self.df) * dr
        else:
            raise ValueError("mode can be bounded or free")

        if inplace:
            self.df = newdf.copy()
        else:
            return newdf

    def add_gaussian_noise(self, sigma=0.1, inplace=False):
        """add gaussian noise to the data. Results may be negative or above 1"""
        # random.normal accepts x and takes its shape to return random values.
        # so, the addition of x and the ouptut of np.random.normal is as
        # expected: points by point.
        newdf = self.df.apply(lambda x:x+np.random.normal(x, scale=sigma))
        if inplace:
            self.df = newdf.copy()
        else:
            return newdf

    def _make_df_compatible(self):
        """if you use MakeBuilder, you can really add any combi of experiments
        However, the resulting df built is not MIDAS compatible for sure. For example,
        not all same time are available for each experiment.

        """
        raise NotImplementedError

    def get_residual_errors(self, level="time", normed=False):
        """Return vector with residuals errors

        The residual errors are interesting to look at in the context
        of a boolean analysis. Indeed, residual errors is the minimum error
        that is unavoidable with a boolean network and comes from the discrete
        nature of such a model. In a boolean analysis, one would compare 0/1
        values to continuous values between 0 and 1. Therefore, however good
        is the optimisation, the value of the goodness of fit term cannot go
        under this residual error.

        :param: level to sum over.
        :return: returns residual errors :math:`\sum (round(x)-x)^2`
            the summation is performed over species and experiment by default

        .. doctest::

            >>> from cellnopt.core import cnodata, XMIDAS
            >>> m = XMIDAS(cnodata("MD-ToyMMB_T2.csv"))
            >>> m.get_residual_errors()
            time
            0       0.000000
            10      2.768152
            100     0.954000
            dtype: float64

        """
        #FIXME use normed to divide residual errors by appropriate N
        diff = (self.df - self.df.apply(lambda x: x.round(), axis=1))
        diff_square = diff.apply(lambda x: x**2, axis=1)
        S = diff_square.sum(level="time").sum(axis=1)

        # FIXME. do we take time 0 into account ?
        if normed :
            S /= len(self.experiments) * len(self.times)
        return S

    def copy(self):
        x = XMIDAS()
        x._data = self._data.copy()
        # FIXME
        #x._missing_time_zero = self._missing_time_zero
        x.cellLine = self.cellLine
        x.df = self.df.copy()
        x._experiments = self.experiments.copy()
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

    def correlation_experiment_one_signal(self, name):
        t = self.df[name]
        pylab.pcolor(t.unstack(1).corr())

    def pca(self, signal, pca_components=2):
        """Not sure this is the proper way...

        get all experiment related to 1 signal

        .. plot::
            :include-source:
            :width: 80%

            from cellnopt.core import *
            m = midas.XMIDAS(cnodata("MD-ToyPB.csv"))
            #m.df = abs(m.df)
            #m.df/=m.df.max()
            m.pca("gsk3")


        .. todo::  pls = PLSRegression(n_components=3)
            from sklearn.pls import PLSCanonical, PLSRegression

        """
        from sklearn.decomposition import PCA
        pca = PCA(n_components=pca_components)
        #t = self.df[signal]
        #pca.fit(t.unstack(1))

        t = self.df[signal].fillna(0)
        X = t.unstack()
        X_r = pca.fit(X).transform(X)
        pylab.plot(X_r[:,0], X_r[:,1], 'o')
        print(pca.explained_variance_ratio_)

        return pca

    def boxplot(self, mode="time"):
        """

        :param str mode: time or species

        .. plot::

            from cellnopt.core import *
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.boxplot(mode="other")
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

            from cellnopt.core import *
            m = XMIDAS(cnodata("MD-ToyPB.csv"))
            m.radviz(["ap1", "gsk3", "p38"])

        """
        if species == None:
            species = list(self.df.columns)
        from pandas.tools.plotting import radviz
        df = self.df.reset_index()
        del df['time']
        del df['cellLine']
        pylab.figure(1)
        pylab.clf()
        radviz(df[['experiment']+species], "experiment")
        pylab.legend(fontsize=fontsize)

    def discretize(self, **kargs):
        return self.discretise(**kargs)

    def discretise(self, inplace=True,N=2):
        """

        :param int N: number of discrete values (defaults to 2). If set to
            2, values will be either 0 or 1. If set to 5, values wil lbe in
            [0,0.25,0.5,0.75,1]
        :param inplace:

        .. warning. data has to be normalised
        """
        assert N>=1
        N = N-1.
        self.logging.info("Discretization between 0-1 assuming normalised data")
        if inplace:
            self.df = self.df.apply(lambda x: (x*N).round()/N)
        else:
            df = self.df.apply(lambda x: (x*N).round()/N)
            #df = df.values.round(1)
            return df

    def round(self, inplace=True, decimals=0):
        if inplace == True:
            self.df = self.df.apply(lambda x : x.round(decimals=decimals))
        else:
            return self.df.values.apply(lambda x : x.round(decimals=decimals))

    def hcluster(self, mode="experiment"):
        """


        .. plot::
            :include-source:
            :width: 80%

            from cellnopt.core import *
            m = midas.XMIDAS(cnodata("MD-ToyPB.csv"))
            m.hcluster("species")

        """
        assert mode in ["experiment", "time", "species"]
        from scipy.spatial.distance import pdist, squareform
        from scipy.cluster.hierarchy import linkage, dendrogram
        pylab.clf()
        if mode == "experiment":
            distxy = squareform(pdist(self.df.unstack("time"), metric='euclidean'))
        elif mode == "time":
            distxy = squareform(pdist(self.df.unstack("experiment"), metric='euclidean'))
        elif mode == "species":
            distxy = squareform(pdist(self.df.transpose(), metric='euclidean'))
        R = dendrogram(linkage(distxy, method='complete'))


        if mode == "time":
            pylab.xticks(pylab.xticks()[0], self.times)
            pylab.title("Clustering by time")
        elif mode == "experiment":
            pylab.xticks(pylab.xticks()[0], self.experiments.index)
            pylab.title("Clustering by experiments")

        else:
            pylab.xticks(pylab.xticks()[0],self.species)
            pylab.title("Clustering by species")

    def heatmap(self, cmap="heat", transpose=False):
        """Hierarchical clustering on species and one of experiment/time level

        .. plot::
            :include-source:
            :width: 80%

            from cellnopt.core import *
            m = midas.XMIDAS(cnodata("MD-ToyPB.csv"))
            m.heatmap()

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

        :param str mode: `timeseries` shuffles experiments and species; timeseries
            are unchanged. `all` shuflles through time, experiment and species.


        mode can be

        # `timeseries` that is
        # `all`
        # `signals` or `species`: sum over signals is constant


        #. by_signals (or by_species, by_columns, species, signals,
            columns) shuffles each column independently. All values are shuffled
            but the sum over a column/species remains identical.
            constqnt is df.sum()
        # shuffle over index. This means that values with same cell/exp/time are shuffled;
          This is therefore over species as well but keep a kind of time information
          constqnt is sum over experiment: m.df.sum(level="experiment").sum(axis=1)


        .. plot::
            :width: 80%

            from cellnopt.core import *
            m = midas.XMIDAS(cnodata("MD-ToyPB.csv"))
            m.plot()

        Shuffling qll timeseries keeping their structures:

        .. plot::
            :include-source:
            :width: 80%

            from cellnopt.core import *
            m = midas.XMIDAS(cnodata("MD-ToyPB.csv"))
            m.shuffle(mode="timeseries")
            m.plot()


        """
        if inplace != True:
            raise NotImplementedError

        if mode == "experiment":
            self.df.reindex(self.df.index, key=lambda x:
                list(self.experiments.index).index(x[1]))
        elif mode == "all":
            # The random.shuffle function does not work!! somehow sum of data increases or decreases
            # One must use numpy.random.shuffle instead
            #print(shuffle)
            shape = self.df.shape
            data = self.df.values.reshape(shape[0]*shape[1])
            np.random.shuffle(data)
            count = 0
            # not very efficient but works for now
            for i in range(0,shape[0]):
                for j in range(0,shape[1]):
                    self.df.values[i][j] = data[count]
                    count += 1
        elif mode in ["signals", "species",  "columns"]:
            for c in self.df.columns:
                self.df[c] = np.random.permutation(self.df[c].values)
        elif mode == "indices":
            # m.df.sum(level="experiment").sum(axis=1) is constant
            for this_index in self.df.index:
                np.random.shuffle(self.df.ix[this_index].values)
        elif mode == "timeseries":
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
            raise NotImplementedError

    # works for simple cases where only one stimuli is on at a time
    def sort_experiments_by_stimuli(self):
        list_exps = []

        for stimulus in self.names_stimuli:
            list_exps.append(('Stimuli', stimulus))
        for inhibitor in self.names_inhibitors:
            list_exps.append(('Inhibitors', inhibitor))
        new_order = self.experiments.groupby(list_exps).groups.values()
        new_order = list(pylab.flatten(new_order))
        print(new_order)
        self._experiments = self.experiments.reindex_axis(new_order, axis=0)

    def sort_experiments_by_inhibitors(self):
        list_exps = []
        for inhibitor in self.names_inhibitors:
            list_exps.append(('Inhibitors', inhibitor))
        for stimulus in self.names_stimuli:
            list_exps.append(('Stimuli', stimulus))

        new_order = self.experiments.groupby(list_exps).groups.values()
        new_order = list(pylab.flatten(new_order))
        print(new_order)

        self._experiments = self.experiments.reindex_axis(new_order, axis=0)

    def __eq__(self, other):
        if all(other.df == self.df) == False:
            return False
        if all(other.experiments == self.experiments) == False:
            return False
        return True

"""
class TrendTimeSeries(object):
    def __init__(self, data=None, times=None):
        if isinstance(data, pd.TimeSeries):
            self.data = data

        # pd.DataFrame([0,2,4,6,8], index=[0,10,20,30,40])

    def _set_times(self):
        pass
    def _get_times(self):
        pass
    times = property(_set_times, _get_times, doc="")

    def _set_data(self):
        # if a timeseries ok
        # otherwise transform to timeseries
        pass
    def _get_data(self):

        pass
    data = property(_get_data,_set_data)
"""


class TypicalTimeSeries(object):
    """Utility that figures out the trend of a time series

    Returns color similar to what is contained in DataRail.


    .. todo:: must deal with NA


    """
    def __init__(self, times=None):
        self._times = times  # ref do not change
    def _get_times(self):
        return self._times
    times = property(_get_times)

    def transient(self, x=None):
        """

        m = MIDASReader(...)
        y = transient(m.times)
        x = m.times
        plot(x,y)

        returns normqlised vector
        """
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

    def plot(self, data):
        corrs = self._get_correlation(data)
        clf()
        pylab.plot(self.times, self._normed(data), label="data", lw=2, ls="--")
        # transient
        pylab.plot(self.times, self.transient(), 'o-', label="transient " + str(corrs['transient']))
        # earlier
        pylab.plot(self.times, self.earlier(), 'o-', label="earlier " + str(corrs['earlier']))
        pylab.plot(self.times, self.earlier(n=1, N=10), 'o-', label="earlier2 " + str(corrs['earlier2']))
        # later
        pylab.plot(self.times, self.later(), 'o-', label="later " + str(corrs['later']))
        # constant
        pylab.plot(self.times, self.constant(.5), 'o-', label="constant " + str(corrs['constant_half']))
        # sustained
        pylab.plot(self.times, self.sustained(L=.5), 'o-', label="sustained" + str(corrs['sustained']))
        pylab.plot(self.times, self.inverse_sustained(L=.5), 'o-', label="inv sustained" + str(corrs['inverse_sustained']))
        pylab.legend()

    def get_bestfit(self, data):
        corrs = self._get_correlation(data)
        keys,values = (corrs.keys(), corrs.values())
        M  = max(values)
        return keys[np.argmax(values)]

    def get_bestfit_color(self, data):
        corrs = self._get_correlation(data)
        keys,values = (corrs.keys(), corrs.values())
        M  = max(values)
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


