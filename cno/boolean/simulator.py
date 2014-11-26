import random

import numpy
import pylab
import pandas as pd

__all__ = ["RandomSimulator", "Simulator"]


# Set a different ordering to deal with nan values
def _bool3key(x):
    """
    Defines the keys used to order the list.
    The only allowed values are True, False, 1,0 and numpy.nan.
    """ 
    return _bool3key.__logic_sort__[x]
_bool3key.__logic_sort__ = {0:-1,numpy.nan:0,1:1}

def and3(*args):
    return min(*args,key=_bool3key)

def or3(*args):
    return max(*args,key=_bool3key)

    
class Simulator(object):
    def __init__(self):
        self._values = pd.DataFrame()

    def _get_values(self):
        return self._values
    def _set_values(self, df):
        self._values = df
        return self._values
    values = property(_get_values, _set_values, doc="")

    def plot_time_course(self, data, mode='boolean', fontsize=16):
        # TODO sort columnsi alphabetically
        # FIXME: twiny labels are slightly shifted
        # TODO flip 
        if mode == 'boolean':
            cm = pylab.get_cmap('gray')
            pylab.clf() 
            data = pd.DataFrame(data).fillna(0.5)
            pylab.pcolor(data, cmap=cm, vmin=0, vmax=1, 
                shading='faceted')
            pylab.colorbar()
            ax1 = pylab.gca()
            ax1.set_xticks([])
            Ndata = len(data.columns)
            ax1.set_xlim(0, Ndata)
            ax = pylab.twiny()
            ax.set_xticks(pylab.linspace(0.5, Ndata+0.5, Ndata ))
            ax.set_xticklabels(data.columns, fontsize=fontsize, rotation=90)
            times = list(data.index)
            Ntimes = len(times)
            ax1.set_yticks([x+0.5 for x in times])
            ax1.set_yticklabels(times[::-1],  
                    fontsize=fontsize)
            pylab.sca(ax1)
        else:
            print('not implemented')



class _Simulator(object):
    """A base class for diverse simulators based on boolean logic"""

    def __init__(self, cnograph):
        self.cnograph = cnograph
        self.cues = self.cnograph.midas.experiments
        
        self.init()

        #self.stimuli = list(self.cnograph.midas.stimuli.columns)
        #self.inhibitors = list(self.cnograph.midas.inhibitors.columns)
        #self.inhibitors = [x[:-1] for x in self.inhibitors]
        #self.nameCues = self.stimuli + self.inhibitors

        

    def _get_stimuli(self):
        return list(self.cnograph.midas.stimuli.columns)
    stimuli = property(_get_stimuli)

    def _get_inhibitors(self):
        inhibitors = list(self.cnograph.midas.inhibitors.columns)
        return [x[:-1] for x in inhibitors]
    inhibitors = property(_get_inhibitors)
    
    def _get_cues(self):
        return self.stimuli + self.inhibitors
    nameCues = property(_get_cues)
    
    def init(self):
        """create some useful attributes to store

	* node names
	* edges
        * predecessors
        * data

	"""
        self.tick = 0
        self.predecessors = {}
        self.nodes = self.cnograph.nodes()[:]

        for node in self.cnograph.nodes():
            self.predecessors[node] = self.cnograph.predecessors(node)

        self.edges = {}
        for k in self.cnograph.edge.keys():
            self.edges[k] = self.cnograph.edge[k]


        self.data = {} # data on each node
        self.dataplus = {}
        for k in self.nodes:
            self.data[k] = [numpy.nan]
            self.data[k] = [0]

        for edge in self.edges.keys():
            for node in self.edges[edge]:
                self.edges[edge][node]["weight"] = 1

        self.links = {}
        for node in self.nodes:
            self.links[node] = [self.edges[pred][node]['link'] for pred in self.predecessors[node]]
        self.weights = {}
        for node in self.nodes:
            self.weights[node] = [self.edges[pred][node]['weight'] for pred in self.predecessors[node]]
            
        self.species = sorted(self.normalNodes())

    def step_simulation(self):
        raise NotImplementedError

    def set_initial_values(self, values={}):
        for k,v in values.items():
            self.data[k][self.tick] = v
 
    def simulate(self, N=10, reset=True):
        if reset: self.init()
        for i in range(0,N):
            self.step_simulation()

    def preprocessing(self):
        self.cnograph.preprocessing()
        self.init()

    def andNodes(self):
        return [x for x in self.nodes if "^" in x]

    def normalNodes(self):
        return [x for x in self.nodes if "^" not in x and x not in self.stimuli]

    def plot(self):
        """
        
        .. plot::
            :include-source:
            :width: 80%
            
            from cellnopt.simulate import *
            from cellnopt.core import *
            pkn = cnodata("PKN-ToyPB.sif")
            midas = cnodata("MD-ToyPB.csv")
            s = boolean.BooleanSimulator(CNOGraph(pkn, midas))
            s.simulate(30)
            s.plot()
        """
        
        pylab.clf()

        data = numpy.array([self.data[x] for x in self.species if x in self.species])
        data = data.transpose()
        data = 1 - pylab.flipud(data)

        pylab.pcolor(data, vmin=0, vmax=1, edgecolors="k")
        pylab.xlabel("species"); 
        pylab.ylabel("Time (tick)");
        pylab.gray()

        pylab.xticks([0.5+x for x in range(0,30)], self.species, rotation=90)

        pylab.ylim([0, self.tick])
        pylab.xlim([0, len(self.species)])

    def plotTimeCourses(self, species, hold=True):
        if isinstance(species, str) == True:
            species = [species]
        if hold==False:
            pylab.clf()

        for this in species:
            x = self.data[this]
            N = len(x)
            pylab.plot(x, "o-", label=this); 
            pylab.hold(True)
        pylab.legend()
        pylab.axis([0, N, -0.1, 1.1]); 
        pylab.grid()




class RandomSimulator(Simulator):
    def __init__(self, cnograph):
        super(RandomSimulator, self).__init__(cnograph)

    def step_simulation(self):
        for this in self.normalNodes():
            value = random.random()
            self.data[this].append(value)
        self.tick += 1


#def test(initvalues=None):
#    s = Simulator(CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv")))
#    s.stimuli = ["tnfa", 'egf']
#    s.set_initial_values({"egf":1, "tnfa":1})
#    s.simulate(30)
#    s.plot()
#    return s




