# -*- python -*-
#
#  This file is part of the cinapps.tcell package
#
#  Copyright (c) 2012-2013 - EMBL-EBI
#
#  File author(s): Thomas Cokelaer (cokelaer@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  cellnopt.core website: http://www.cellnopt.org
#
##############################################################################
"""Module dedicates to the metabolites CNA format

:Status: mature but not all features implemented.

Based on load_substances_inter matlab version provided by Steffen, Klamt, MPI 
"Dynamik komplexer technischer Systeme" Magdeburg

"""
__all__ = ["Metabolites"]


class Metabolites(object):
    """Read metabolites file and convert to a Metabolites data structure.

    Metabolites format is a CSV format that looks like::

        abl         abl         NaN     0   188  380    1    1
        akap79      akap79      NaN     0   989  442    1    1


    Columns are

    #. identifier of this substance in CNA (e.g.: G6P) without blanks 
       stored in :ref:`species`
    #. the full name of the metabolite (e.g.: glucose-6-phosphate; no blanks
        allowed!) stored in :attr:`specLongNames`
    #. could be # or a value stored in :ref:`specDefault`
    #. a value 0/1  (ignored right now)
    #. 4 following columns must be numerical values stored in :attr:`specBoxes`

    .. todo:: specify precicely the content of the columns

    ::

        >>> from cno.io import Metabolites
        >>> m = Metabolites("metabolites")
        >>> m.species


    """
    def __init__(self, filename, verbose=True):
        """.. rubric:: Constructor


        :param str filename: input filename containing the metabolites data
        :param bool verbose: (True by default)


        """

        self.filename = filename
        self.verbose = verbose
        #todo: check that the file exist

        #: attribute to store the species
        self.species = []    # todo: transform to dictionary
        #: attribute to store the specLongNames
        self.specLongNames = []
        #: attributes to store specNotes
        self.specNotes = {}
        #: attribute to store spec default
        self.specDefault = []
        #: attribute to store specBoxes
        self.specBoxes = []
        #ws = char([9:13 32]); %A# whitespace for octave compatibility

        self._read_data()

    def _read_data(self):
        if self.verbose:
            print(' ')
            print('Reading Species ...')

        # scanning the entire file
        f = open(self.filename,'r')
        data = []                   # the data structure to populate
        for line in f.readlines():  # for each line
            # convert tab.to white space, remove trailing and \n character
            line = line.replace('\t',' ').replace('\n','').strip()

            # do not consider commented or empty lines
            if line.startswith("%") or line.startswith('#'): 
                continue
            if len(line) == 0:
                print("Found an empty line. Skipped")
            else:
                data.append(line)
        f.close()

        # The actual processing is done here
        for i, x in enumerate(data):
            # split the line using white space delimiter
            x = x.split()
            self.species.append(x[0])        # Store the species (first column)
            self.specLongNames.append(x[1]) # Store the second column (long name)

            # Then, it depends....
            if x[2] == '#':
                self.specDefault.append('NaN')
            else:
                try:
                    self.specDefault.append(float(x[2]))
                except:
                    self.specDefault.append(x[2])

            # x[3] is skipped... do not know why


            try:
                xpos, ypos, map_nr, rtype = x[4:]
                self.specBoxes.append([i+1, float(xpos), float(ypos), 0, float(map_nr), float(rtype)])
            except:
                self.specBoxes.append([i+1, 50,50,0,1,1])
                print('warning. set default values')
                #raise ValueError('Could not parse line. Missing data/column ?')


        # Some sanity checks.
        check = set(self.species)
        if len(check) != len(self.species):
            raise ValueError('Found a duplicated metabolite name ! Fix the input file.')

        self.N = len(self.species)
        print("Found %s species" % len(self)) 
        # check that all attributes are correct.
        # in principle there is a metabolite_notes file to be read. If not
        # found, should raise a warning.
        self._set_notes()

    def __len__(self):
        return len(self.species)

    def _set_notes(self, filename=None):
        """If a file called metabolites_notes or filename_notes.

        If not found, the attribute specNotes is kept empty
        """
        if filename == None:
            filename = self.filename + "_notes"
        try:
            f = open(filename, "r")
        except IOError:
            # no valid file
            print("No valid notes were found. Skipping.")
            return
        except Exception:
            raise Exception

        # if the file was found, read it and populated specNotes
        [self.specNotes.setdefault(x) for x in self.species]
        for line in f.readlines():
            spec = line.split()[0].strip()
            note = " ".join(line.split()[1:])
            self.specNotes[spec] = note
        f.close()

        assert sorted(self.species) == sorted(self.specNotes.keys())

