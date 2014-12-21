import sys
from cno.apps.plotmodel import plotmodel



# todo
class _test_something(object):
    def setUp(self):
        sys.argv[1] = '--model'
        del sys.argv[2] # remember that -s is in sys.argv[2], see below
        self.args = sys.args
    def test_method(self):
        plotmodel(self.args)


