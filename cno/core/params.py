import argparse
from easydev import AttrDict
import functools

__all__ = ['Parameter', 'Parameters', 'ParamsGA', 'params_to_update',
           'ParamsPreprocessing']


# Let us create a handy named tuple to define a parameter data structure
from collections import namedtuple
Params = namedtuple('Params',
        ['name', 'argname', 'default', 'description'])


def params_to_update():
    """
    Decorator that provides the wrapped function with an attribute 'actual_kwargs'
    containing just those keyword arguments actually passed in to the function.
    """
    def decorator(function):
        @functools.wraps(function)
        def inner(self, *args, **kwargs):
            inner.actual_kwargs = kwargs
            return function(self, *args, **kwargs)
        return inner
    return decorator


class Parameter(Params):
    """Define a user parameter

    A parameter is defined by
    - **name** a **destination** name, which is the variable name to
      be used internally.
    - **argname** a user name to be used within the configuration file. It can
      be used also within executable. It may be longer than **name** to be
      more explicit (with dashes, which are not accepted in Python)
    - **default** a default value
    - **description** to be used with a --help option
    - **types** may be provided.

        >>> p = Parameter('optimisation', '--max-iteration', 100,
            'maximum number of iterations')
        >>> p.value
        100
        >>> p.value = 20
        >>> p.value
        20
        >>> p.reset()
        >>> p.value
        100
        >>> p1 == p2

    Here, we can further check types
    Tuple so nothing can be modified except new arguments such as types
    and value.

    """
    def __init__(self, name, argname, default, description, types=[]):
        self._value = None
        assert argname.startswith("--"), "argument name must start with -- signs"
        super(Parameter, self).__init__(name=name, argname=argname,
                default=default, description=description)
        self.value = default
        # TODO check types

    def _get_value(self):
        return self._value
    def _set_value(self, value):
        self._value = value
    value = property(_get_value, _set_value)

    def reset(self):
        self._value = self.default

    def __eq__(self, other):
        if self.value != other.value:
            return False
        for arg in ['name', 'argname', 'default', 'description']:
            if getattr(self, arg) != getattr(other, arg):
                return False
        return True


class Parameters(AttrDict):
    reserved = ['name', 'description']
    def __init__(self, name, description=""):
        super(Parameters, self).__init__()
        # reserved name
        self.description = description
        self.name = name

    def remove_parameter(self, name):
        del self[name]

    def add_parameter(self, params):
        if params.name in self.reserved:
            raise ValueError("The name of a parameter cannot be one of the reserved word %s" % self.reserved)
        self[params.name] = params

    def _get_names(self):
        # + ['value'] is required. see tests
        return [x for x in self.keys() if x not in self.reserved + ['value']]

    def __repr__(self):
        txt = "Section/Parameters: %s\n" % self.name
        for key in self._get_names():
            txt += "- " + key + ": " + self[key].__repr__() + "\n"
        return txt

    def __str__(self):
        txt =""
        #txt = "; " + self.description
        txt += "[" + self.name + "]\n"
        for k in sorted(self._get_names()):
            try:
                value = self[k].value
            except:
                print "------", k, value
            if value is None:
                value = ""
            elif value is True:
                value = 'True'
            elif value is False:
                value = 'False'
            txt += " = ".join([k, str(value)])
            txt += "\n"
        return txt

    def reset(self):
        for name in self._get_names():
            self[name].reset()

    def __eq__(self, other):
        if sorted(self._get_names()) != sorted(other._get_names()):
            print("False 1")
            return False
        for key in other._get_names():
            p1 = other[key]
            p2 = self[key]
            if (p1 == p2) is False:
                print("False 2" + key)
                return False
        return True


class ParamsGeneral(Parameters):
    def __init__(self):
        super(ParamsGeneral, self).__init__('General', "General description")
        self._init()

    def _init(self):
        for data in [
            ("pknmodel", "--pkn-model", None,
            "The input prior knowledge network in SIF format (or SBML-qual)"),
            ("data", "--data", None, "The input data in MIDAS format"),
            ("formalism", "--formalism", None, "not used yet"),
            ("tag", "--tag", None, "not used yet"),
            ("report", "--report", True,  "create report"),
            ("onweb", "--on-web", True, "open report in a browser. This option also set --report "),
            ("verbose", "--verbose", True,  "verbosity"),
            ("verboseR", "--verbose-R", True,  "verbosity of R scripts")]:

            self.add_parameter(Parameter(*data))


class ParamsPreprocessing(Parameters):
    def __init__(self):
        super(ParamsPreprocessing, self).__init__('Preprocessing', "Preprocessing description")
        self._init()

    def _init(self):
        for data in [("cutnonc", "--with-cutnonc", True, ""),
                     ("compression", "--with-compression", True, ""),
                     ("expansion", "--with-expansion", True, ""),
                     ("maxInputsPerGate", "--max-inputs-per-gate", 3, "")]:
            self.add_parameter(Parameter(*data))


class ParamsGA(Parameters):
    error_msg = {
            'elitism': "elitism must be strictly positive and less than "
            "popsize argument"
        }

    def __init__(self):
        super(ParamsGA, self).__init__('GA', 'Genetic algorithm')
        self._init()

    def _init(self):
        pass
        # adding all info required
        self.add_parameter(Parameter('elitism', '--elitism', 5,
            "The elitism number (should be 10%% of the popsize)"))
        self.add_parameter(Parameter('sizefactor', '--size-factor', 0.0001,
            "The penalty factor (if NaN values)"))
        self.add_parameter(Parameter("popsize", "--population-size", 50,
            "The population size"))
        self.add_parameter(Parameter('maxtime', "--max-time", 60,
            "Maximum time of the simulation (seconds)"))
        self.add_parameter(Parameter('NAFac', "--na-factor", 1,
            "The penalty factor (if NaN values)"))
        self.add_parameter(Parameter('pmutation', "--pmutation", 0.5,
            "Mutation rate"))
        self.add_parameter(Parameter("maxgens", "--max-generations", 500,
            "maximum number of generation"))
        self.add_parameter(Parameter('maxstallgens', "--max-stall-generations", 100,
            "Max number of stall generation"))
        self.add_parameter(Parameter('selpress', "--selection-pressure", 1.2,
            "todo"))
        self.add_parameter(Parameter('reltol', "--relative-tolerance", 0.1,
            "todo"))
        self.add_parameter(Parameter('ga_verbose', "--ga-verbose", True,
            "verbosity in genetic algorithm"))
        # indices starts at zero.
        self.add_parameter(Parameter('time_index_1', "--time-index-1", 1,
           "first time index to optimise"))
        self.add_parameter(Parameter('time_index_2', "--time-index-2", 2,
           "second time index to optimise"))




class CNOConfigParser(AttrDict):

    def __init__(self, filename=None):
        super(CNOConfigParser, self).__init__()
        if filename:
            self.read(filename)

    def remove_section(self, section):
        del self[section]

    def add_section(self, section):
        self[section.name] = section

    def read(self, filename):
        # clean the sections
        for section in self.keys():
            del self[section]

        from ConfigParser import ConfigParser
        config = ConfigParser()
        config.read(filename)
        for this in config.sections():
            sec = eval("Params" + this + "()")
            options = config.options(this)
            for option in options:
                value = config.get(this, option)
                setattr(sec, 'value', value)
            self[this] = sec

        return config

    def to_ini(self, filename):
        self.save(filename)

    def save(self, filename):
        with open(filename, 'w') as fh:
            txt = self.__str__()
            fh.write(txt)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        txt = ""
        for k in self.keys():
            section = self[k]
            txt += section.__str__()
            txt += "\n"
        return txt

    def reset(self):
        for section in self.keys():
            self[section].reset()

    def __eq__(self, other):
        if sorted(other.keys()) != sorted(self.keys()):
            print('config False 1')
            return False
        for key in other.keys():
            # check that all sections are identical
            if (other[key] == self[key] ) is False:
                print('config False 2 ' + key)
                return False
        return True


class OptionsBase(argparse.ArgumentParser):
    def __init__(self, version="1.0", prog="cellnopt"):
        super(OptionsBase, self).__init__(version=version, prog=prog)
        from cno.core.params import CNOConfigParser
        self.config = CNOConfigParser()

        section = ParamsGeneral()
        self.add_section(section)

        section = ParamsPreprocessing()
        self.add_section(section)

    def add_section(self, section):
        self.config.add_section(section)
        group = self.add_argument_group(section.name, section.description)
        for key in sorted(section._get_names()):
            param = section[key]

            if isinstance(param.default, bool):
                action = "store_" + str(param.default).lower()
                group.add_argument(param.argname, dest=param.name,
                    action=action, help=param.description)
            elif isinstance(param.default, str):
                group.add_argument(param.argname, dest=param.name,
                    default=param.default, type=str,
                    help=param.description)
            else:
                group.add_argument(param.argname, dest=param.name,
                    default=param.default, help=param.description)


class CNOConfig(CNOConfigParser):

    def __init__(self, filename=None):
        super(CNOConfig, self).__init__(filename)
        self.init_config()

    def init_config(self):
        self.add_section(ParamsGeneral())
        self.add_section(ParamsGA())
        self.add_section(ParamsPreprocessing())







