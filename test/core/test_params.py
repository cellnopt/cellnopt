from cno.core.params import Parameter, Parameters
from cno.core.params import ParamsGA, ParamsGeneral, ParamsPreprocessing
from cno.core.params import CNOConfigParser, CNOConfig
from cno.core.params import OptionsBase


def test_parameter():
    p = Parameter('opt', '--opt', 20, 'description')
    assert p.value == 20
    assert p.default == 20
    p.value = 10
    assert p.value == 10
    assert p.default == 20
    p.reset()
    assert p.value == 20
    assert p.default == 20


def test_parameter_eq():
    p1 = Parameter('opt', '--opt', 20, 'description')
    p2 = Parameter('opt', '--opt', 20, 'description')
    assert p1 == p2
    p1.value = 10
    assert (p1 != p2) is False


def test_parameters():
    p = Parameters('test', 'a test example')
    param = Parameter('opt', '--opt', 20, 'description')
    p.add_parameter(param)
    p['opt'].value = 10
    print(p)
    p  # test __repr__
    p.reset()
    assert p['opt'].value == 20

    param = Parameter('opt', '--opt', 20, 'description')
    p.remove_parameter('opt')

    try:
        param = Parameter('name', '--name', 20, 'description')
        p.add_parameter(param)
        assert False
    except:
        assert True


def test_parameters_eq():
    p1 = ParamsGeneral()
    p2 = ParamsGeneral()
    p3 = ParamsPreprocessing()
    assert p1 == p2
    assert (p1 != p3) is True

def test_config_parser():
    s1 = ParamsGA()
    s2 = ParamsGeneral()
    c1 = CNOConfigParser()
    c1.add_section(s2)
    c1.add_section(s1)

    s1 = ParamsGA()
    s2 = ParamsGeneral()
    c2 = CNOConfigParser()
    c2.add_section(s2)
    c2.add_section(s1)

    assert c1 == c2

    from easydev import TempFile
    fh = TempFile()
    c1.save(fh.name)
    c2 = CNOConfigParser(fh.name)
    fh.delete()
    assert c1 == c2

def test_cno_config():
    c = CNOConfig()
    c.init_config()
    c.remove_section('GA')

def test_decorator():
    from cno.core.params import params_to_update

    class A(object):
        def __init__(self):
            pass
        @params_to_update
        def runme(self, arg1=1, arg2=2):
            self.runme.actual_kwargs
    a = A()
    a.runme()
    a.runme(arg1=2)
    assert "arg1" in a.runme.actual_kwargs


def test_options():
    o = OptionsBase()
