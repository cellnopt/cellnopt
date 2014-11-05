"""This is an experimental module to design an XML 
format for MIDAS data sets.

"""
import xml.etree.ElementTree as ET

from cno.io.measurements import Measurements
from cno.io.measurements import Measurement

__all__ = ["XMLMIDAS"]


class XMLMIDAS(object):
    """Class to read MIDAS in XML format. """
    def __init__(self, filename):
        self.tree = ET.parse(filename)
        self.root = self.tree.getroot()
        self._data = None

    def readXML(self):
        """Read XML MIDAS file and build up a Measurements instance"""
        # scan all measurements
        mb = Measurements()
        for this in self.tree.findall("measurements")[0].findall("measurement"):

            tr = this.find('treatments')
            species_parser = this.find('species')
            time_parser = this.find('time')
            cell_parser = this.find('cell')

            # get time 
            time = float(time_parser.find('value').text)
            time_units = time_parser.find('value').attrib['units']

            # get cell name
            cellname = cell_parser.find('name').text

            # get species information
            species_name = species_parser.find('name').text
            value = species_parser.find('value').text

            # get stimuli and inhibitors
            dict_stimuli = {}
            dict_inhibitors = {}
            if tr.find('stimuli') is not None:
                for stimulus in tr.find('stimuli').findall('stimulus'):
                    name = stimulus.find('name').text
                    value = float(stimulus.find('value').text)
                    dict_stimuli[name] = value

            if tr.find('inhibitors') is not None:
                for inhibitor in tr.find('inhibitors').findall('inhibitor'):
                    name = inhibitor.find('name').text
                    value = float(inhibitor.find('value').text)
                    dict_inhibitors[name] = value

            measure = Measurement(species_name, time, dict_stimuli, 
                    dict_inhibitors, value, 
                    cellLine=cellname, units=time_units)
            mb.append(measure)
        return mb
            
    def _get_midas(self):
        from cno.io.measurements import MIDASBuilder

        # create measurements 
        measurements = self.readXML()

        # add them to a MIDASBuilder
        xm = MIDASBuilder()
        xm.add_measurements(measurements)

        # get the MIDAS file itself
        midas = xm.xmidas
        return midas
    xmidas = property(_get_midas, doc='Get XMIDAS instance from the XML')

