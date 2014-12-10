"""




"""

__all__ = ['Parameter', 'Parameters', 'BooleanParameters']


class Parameter(object):
    """Define a user parameter

    A parameter is defined by 
    - **dest** a **destination** name, which is the variable name to 
      be used internally. 
    - **name** a user name to be used within the configuration file. It can 
      be used also within executable. The difference with the **dest** name is 
      that it may be longer to be more explicit (with dashes, which are not accepted in Python)
    - **section** a section to store the parameter in the configuration file
    - **default** a default value
    - **decription** to be used with a --help option
    - **type** may be provided.

    p = Parameter('optimisation', '--max-iteration', 'maxiter', 
        100, 'maximum number of iterations')

    """
    def __init__(self, section, name, dest, default, types=None, description=""):
        self.section = section
        self.name = name
        self.dest = dest
        self.default = default
        if types is None:
            self.types = type(default)
        self.description = description
        #self.value = None

    def _get_kargs(self):
        kargs = {
            'dest':self.dest,
            'default':self.default,
            'help': self.description.replace("%", "%%"),
            'type': self.type}
        return kargs

    def __str__(self):
        txt = self.name
        return txt


class Parameters(object):
    def __init__(self):
        self.parameters = {}

    def add_parameter(self, section, name, dest, default,
                      description=""):
        # Could add a Parameter instance ?
        p = Parameter(section, name, dest, default,  description)
        self.parameters[name] = p

    def get_keys_from_section(self, section):
        keys = [k for k,v in self.parameters.iteritems() if v.section == section]
        keys = sorted(keys)
        return keys


class GAParameters(Parameters):
    # THe keys used here have the same caps as in the R code.
    name = 'GA'
    params = {
        "sizeFac": 0.0001,
        "NAFac": 1,
        "popSize": 50,
        "pMutation": 0.5,
        "maxTime": 60,
        "maxGens": 500,
        "stallGenMax": 100,
        "selPress": 1.2,
        "elitism": 5,
        "relTol": 0.1,
        "verbose": True,
        "timeIndex1": 2,
        "timeIndex2": 3
        }
    def __init__(self):
        super(GAParameters, self).__init__()
        self._init_gabinary_t1()

    def _init_gabinary_t1(self):
        # just an alias
        default = self.params
        name = self.name
        # adding all info required
        self.add_parameter(name, "--elitism", "elitism", default['elitism'],
                           "The elitism number (should be 10% of the popsize)")
        self.add_parameter(name, "--size-factor", "sizeFac", default["sizeFac"],
                           "The penalty factor (if NaN values)")
        self.add_parameter(name, "--population-size", "popSize", default["popSize"],
                           "The population size")
        self.add_parameter(name, "--max-time", "maxTime", default['maxTime'],
                           "Maximum time of the simulation (seconds)")
        self.add_parameter(name, "--na-factor", "NAFac", default['NAFac'],
                           "The penalty factor (if NaN values)")
        self.add_parameter(name, "--pmutation", "pMutation", default['pMutation'],
                           "Mutation rate")
        self.add_parameter(name, "--max-generation", "maxGens",  default['maxGens'],
                            "maximum number of generation")
        self.add_parameter(name, "--max-stall-generation", "stallGenMax", default['stallGenMax'],
                           "Max number of stall generation")
        self.add_parameter(name, "--selection-pressure", "selPress", default['selPress'],
                           "todo")
        self.add_parameter(name, "--relative-tolerance", "relTol", default['relTol'],
                            "todo")
        self.add_parameter(name, "--ga-verbose", "verbose", True,
                           "verbosity in GA")
        self.add_parameter(name, "--time-index-1", "timeIndex1",  default['timeIndex1'],
                         "first time index to optimise")
        self.add_parameter(name, "--time-index-2", "timeIndex2",  default['timeIndex2'],
                         "second time index to optimise")
    #def __getattr__(self, key):
    #    return self.parameters[key]

