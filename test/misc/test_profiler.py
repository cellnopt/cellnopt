from cno.misc.profiler import do_profile




def test_profiler():

    @do_profile()
    def longrun():

        for i in range(1,10000):
            x = 5.5**2

    
    longrun()
