# -*- python -*-
#
#  This file is part of CNO software
#
#  Copyright (c) 2013-2014 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: http://github.com/cellnopt/cellnopt
#
##############################################################################
import argparse
from easydev import AttrDict
import functools


__all__ = ['Parameter', 'Parameters', 'CNOConfig', 'ParamsGA', 'ParamsGA2',
           'params_to_update',
           'ParamsPreprocessing', 'ParamsDT', 'ParamsFuzzy', 'ParamsSSM']


# Let us create a handy named tuple to define a parameter data structure
from collections import namedtuple
Params = namedtuple('Params',
        ['name', 'argname', 'default', 'description'])


import wrapt
@wrapt.decorator
def params_to_update(wrapped, instance, args, kwargs):
    """
    Decorator that provides the wrapped function with an attribute 'actual_kwargs'
    containing just those keyword arguments actually passed in to the function.
    """
    vars(wrapped)['actual_kwargs'] = kwargs
    return wrapped(*args, **kwargs)



"""def params_to_update():
    def decorator(function):
        @functools.wraps(function)
        def inner(self, *args, **kwargs):
            inner.actual_kwargs = kwargs
            return function(self, *args, **kwargs)
        return inner
    return decorator
"""

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

    def as_dict(self):
        d =  dict([(k, self[k].value)
            for k in self._get_names()])
        return d

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
                print("False 2 " + key)
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
            ("verboseR", "--verbose-R", True,  "verbosity of R scripts"),
            ("cnodata", "--use-cnodata", True,  "search for data/model in cnodata repository"),
            ("config_file", "--config-file", None,  "todo")
        ]:
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


class ParamsGA2(ParamsGA):
    def __init__(self):
        super(ParamsGA, self).__init__('GA2', 'Genetic algorithm')
        self._init()

    def _init(self):
        # adding all info required
        self.add_parameter(Parameter('elitism', '--elitism2', 5,
            "The elitism number (should be 10%% of the popsize)"))
        self.add_parameter(Parameter('sizefactor', '--size-factor2', 0.0001,
            "The penalty factor (if NaN values)"))
        self.add_parameter(Parameter("popsize", "--population-size2", 50,
            "The population size"))
        self.add_parameter(Parameter('maxtime', "--max-time2", 60,
            "Maximum time of the simulation (seconds)"))
        self.add_parameter(Parameter('NAFac', "--na-factor2", 1,
            "The penalty factor (if NaN values)"))
        self.add_parameter(Parameter('pmutation', "--pmutation2", 0.5,
            "Mutation rate"))
        self.add_parameter(Parameter("maxgens", "--max-generations2", 500,
            "maximum number of generation"))
        self.add_parameter(Parameter('maxstallgens', "--max-stall-generations2", 100,
            "Max number of stall generation"))
        self.add_parameter(Parameter('selpress', "--selection-pressure2", 1.2,
            "todo"))
        self.add_parameter(Parameter('reltol', "--relative-tolerance2", 0.1,
            "todo"))
        self.add_parameter(Parameter('ga_verbose', "--ga-verbose2", True,
            "verbosity in genetic algorithm"))
        # indices starts at zero.
        self.add_parameter(Parameter('time_index_2', "--time-index-2", -1,
           "first time index to optimise"))


class ParamsDT(Parameters):

    error_msg = {}

    def __init__(self):
        super(ParamsDT, self).__init__('DiscreteTime', 'description discrete time')
        self._init()

    def _init(self):
        self.add_parameter(Parameter('bool_updates', "--bool-updates", 10,
           "description todo"))
        self.add_parameter(Parameter('lower_bound', "--lower-bound", 0.8,
           "description todo"))
        self.add_parameter(Parameter('upper_bound', "--upper-bound", 10,
           "description todo"))


class ParamsFuzzy(Parameters):

    error_msg = {}

    def __init__(self):
        super(ParamsFuzzy, self).__init__('Fuzzy', 'description to be done')
        self._init()

    def _init(self):

        self.add_parameter(Parameter('do_refinement', '--do-refinement', True,
                                     'description to do'))
        self.add_parameter(Parameter('optimisation_algorithm', '--optimisation-algorithm',
                                     'NLOPT_LN_SBPLX',
                                     'description to do'))

        self.add_parameter(Parameter('optimisation_xtol_abs', '--optimisation-xtol-abs', 0.001,
                                     'description to do'))
        self.add_parameter(Parameter('optimisation_max_eval', '--optimisation-max-eval', 10000,
                                     'description to do'))
        self.add_parameter(Parameter('optimisation_max_time', '--optimisation-max-time', 300,
                                     'description to do'))
        self.add_parameter(Parameter('N', '--multiple-runs', 2,
                                     'description to do'))


class ParamsSSM(Parameters):

    error_msg = {}

    def __init__(self):
        super(ParamsSSM, self).__init__('SSM', 'description to be done')
        self._init()

    def _init(self):
        self.add_parameter(Parameter("maxtime", "--max-time", 60,
                           "maximum time for the optimisation"))
        self.add_parameter(Parameter("dim_ref_set", "--dim-ref-set", 10,
                           "ssm parameter"))
        self.add_parameter(Parameter('n_diverse', "--n-diverse",10,
                           "ssm parameter"))
        # not to be used
        #self.add_parameter(Parameter("verbose", "--ssm-verbose", True,
        #                   "todo"))
        self.add_parameter(Parameter("transfer_function", "--transfer-function", 3,
                           "number of transfer function to be used by the ODE logical formalism"))
        self.add_parameter(Parameter("reltol", "--relative-tolerance", 1e-4,
                           "todo"))
        self.add_parameter(Parameter("atol", "--absolute-tolerance", 1e-3,
                           "todo"))
        # Inf string is to be kept with this spelling and caps as in R syntax
        self.add_parameter(Parameter("maxeval", "--max-evaluation", 'Inf',
                           "maximum number of evaluation during the optimisation"))
        self.add_parameter(Parameter("maxstepsize", "--max-step-size", 'Inf',
                           "todo optimisation"))

        # TODO
        """
        ode_parameters = NULL, indices = NULL,
            local_solver = NULL,
             time = 1,
             mmaxNumSteps = 1e+05, maxErrTestsFails = 50, \
            nan_fac =1,  useVariances=F,initial_state=0.1)
         """


class CNOConfigParser(AttrDict):

    def __init__(self, filename=None):
        super(CNOConfigParser, self).__init__()
        if filename:
            self.read(filename)

    def remove_section(self, section):
        del self[section]

    def add_section(self, section):
        self[section.name] = section

    def as_dict(self):
        d = {}
        for name in self.keys():
            d[name] = self[name].as_dict().copy()
        return d

    def read(self, filename):
        # clean the sections
        for section in self.keys():
            del self[section]

        from ConfigParser import ConfigParser
        config = ConfigParser()
        config.optionxform = str
        config.read(filename)
        for this in config.sections():
            sec = eval("Params" + this + "()")
            options = config.options(this)
            for option in options:
                value = config.get(this, option)
                # HERE, we need to interpret the content of the config file.
                # values could be int or float, in which case we want to cast the string
                # to a float. If the string is True/False, we want also a cast
                try:
                    try:
                        value = int(value)
                    except:
                        value = float(value)
                except:
                    if isinstance(value, str):
                        if value == 'True':
                            value = True
                        elif value == 'False':
                            value = False
                        elif value.strip() == '':
                            value = None
                setattr(getattr(sec, option), 'value', value)
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
    """An ArgumentParser with a CNOConfig file as attribute"""
    def __init__(self, version="1.0", prog="cellnopt"):
        super(OptionsBase, self).__init__(version=version, prog=prog)

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
    """A minimalist configuration class for the various formalisms

    This includes the :class:`ParamsGeneral <General>` section,
    the Genetic Algorithm section and
    the preprocessing section.

    """
    def __init__(self, filename=None):
        super(CNOConfig, self).__init__(filename)
        self.init_config()

    def init_config(self):
        self.add_section(ParamsGeneral())
        self.add_section(ParamsGA())
        self.add_section(ParamsPreprocessing())
