from cno.io.net import NET
from cno import getdata



class testNET(object):

    @classmethod
    def setup_class(klass):
        # test the constrcutors
        klass.net1 = NET()
        klass.net2 = NET(getdata("PKN-ToyMMB.net"))

        try:
            n = NET("dummy")
            assert False
        except:
            assert True

    def test_add(self):
        self.net1.add_net("a -> b +")
        self.net1.add_net("a -> c -")

    def test_accessors(self):
        self.net2.to_sif()
        self.net2.reactions

    def test_print(self):
        print(self.net2)

    def test_write2sif(self):
        import tempfile
        f = tempfile.NamedTemporaryFile()
        self.net2.to_sif(f.name)
        f.close()

    def test_wrong_add_netnet(self):
        try:
            self.net1.add_net("dummy")
            assert False
        except:
            assert True

        try:
            self.net1.add_net("A -> B+")
            assert False
        except:
            assert True

    def test_net2reaction(self):
        assert self.net1.net2reaction("A -> B +") == "A=B"
        assert self.net1.net2reaction("A -> B -") == "!A=B"


