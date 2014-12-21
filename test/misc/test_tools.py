from cno.misc.tools import CNOError



def test_cnoerror():

    try:
        raise CNOError('txt')
        assert False
    except CNOError:
        assert True
