import easydev
import json
# could be prt of easydev itself but for now, let us keep it here
# inside bioservices

__all__ = ['DevTools']


class DevTools(object):
    """wrapper around useful functions.

    See easydev documentation for details.
    """
    def check_range(self, a, b, value):
        easydev.check_range(a, b, value, strict=False)

    def check_param_in_list(self, param, valid_values):
        """

        transforms valid_values into list (e.g., convenient if
        an iterator)
        """
        param = self.tolist(param)
        for name in param:
            easydev.check_param_in_list(name, list(valid_values))

    def swapdict(self, d):
        return easydev.swapdict(d)

    def tolist(self, query):
        """
        'a' ->['a']
        1 -> [1]
        """
        return easydev.codecs.tolist(query)

    def list2string(self, query, sep=",", space=False):
        return easydev.codecs.list2string(query, sep=sep, space=space)

    def to_json(self, dictionary):
        return json.dumps(dictionary)


