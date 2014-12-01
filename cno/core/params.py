

__all__ = ['Parameter', 'Parameters']


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


