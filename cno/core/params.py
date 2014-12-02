

__all__ = ['Parameter', 'Parameters', 'BooleanParameters']


class Parameter(object):
    def __init__(self, section, name, dest, default, description="TODO"):
        self.description = description
        self.section = section
        self.name = name
        self.default = default
        self.dest = dest
        self.type = type(default)
        self.value = None

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
                      description="TODO"):
        p = Parameter(section, name, dest, default,  description)

        self.parameters[name] = p

    def get_keys_from_section(self, section):
        keys = [k for k,v in self.parameters.iteritems() if v.section == section]
        keys = sorted(keys)
        return keys


class BooleanParameters(Parameters):
    # THe keys used here have the same caps as in the R code.
    gaBinaryT1_params = {
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
        super(BooleanParameters, self).__init__()
        self.init_gabinary_t1()

    def init_gabinary_t1(self):
        # just an alias
        default = self.gaBinaryT1_params

        # adding all info required
        self.add_parameter("GA", "--elitism", "elitism", default['elitism'],
                           "The elitism number (should be 10% of the popsize)")
        self.add_parameter("GA", "--size-factor", "sizeFac", default["sizeFac"],
                           "The penalty factor (if NaN values)")
        self.add_parameter("GA", "--population-size", "popSize", default["popSize"],
                           "The population size")
        self.add_parameter("GA", "--max-time", "maxTime", default['maxTime'],
                           "Maximum time of the simulation (seconds)")
        self.add_parameter("GA", "--na-factor", "NAFac", default['NAFac'],
                           "The penalty factor (if NaN values)")
        self.add_parameter("GA", "--pmutation", "pMutation", default['pMutation'],
                           "Mutation rate")
        self.add_parameter("GA", "--max-generation", "maxGens",  default['maxGens'],
                            "maximum number of generation")
        self.add_parameter("GA", "--max-stall-generation", "stallGenMax", default['stallGenMax'],
                           "Max number of stall generation")
        self.add_parameter("GA", "--selection-pressure", "selPress", default['selPress'],
                           "todo")
        self.add_parameter("GA", "--relative-tolerance", "relTol", default['relTol'],
                            "todo")
        self.add_parameter("GA", "--ga-verbose", "verbose", True,
                           "verbosity in GA")
        self.add_parameter("GA", "--time-index-1", "timeIndex1",  default['timeIndex1'],
                         "first time index to optimise")
        self.add_parameter("GA", "--time-index-2", "timeIndex2",  default['timeIndex2'],
                         "second time index to optimise")
    #def __getattr__(self, key):
    #    return self.parameters[key]

