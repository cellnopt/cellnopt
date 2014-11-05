from cno.io.xmlmidas import XMLMIDAS
from cno import getdata

def test_xmlmidas():

    m = XMLMIDAS(getdata("MD-test.xml"))
    xm = m.xmidas
    assert xm.names_stimuli == ['EGF']
    assert xm.names_inhibitors == ['PI3K']



