from cno.misc.models import Models
import pandas as pd
import numpy as np

from easydev import TempFile





def test_models():
    data = np.array([[1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0],
       [1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0],
       [1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0],
       [1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1],
       [1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1],
       [1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0],
       [1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0],
       [1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0],
       [1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0],
       [1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1],
       [1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0],
       [1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1],
       [1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1],
       [1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1],
       [1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1]])

    columns = [u'EGF=PI3K', u'TNFa=PI3K', u'Jnk=cJun', u'PI3K=Akt', u'Raf=Mek', u'!Akt=Mek', u'Mek=p90RSK', u'Mek=Erk', u'Erk=Hsp27', u'TNFa=Jnk', u'TNFa=NFkB', u'TNFa=Hsp27', u'EGF=Raf', u'EGF^TNFa=PI3K', u'Raf^!Akt=Mek', u'Erk^TNFa=Hsp27']

    df = pd.DataFrame(data,columns=columns)
    fh = TempFile()
    df.to_csv(fh.name)
    
    m1 = Models(df)
    m2 = Models(m1)
    m3 = Models(fh.name, index_col=0)  # there is an index column with no name
    fh.delete()

    # trying a stupid constructor
    try:
        Models(1)
        assert False
    except:
        assert True

    assert m1 == m2
    assert m1 == m3

    # plots
    m1.plot()
    m1.plot(1)
    m1.plot('cv')
    m1.errorbar()
    m1.heatmap()

    # exporters
    fh = TempFile()
    m1.to_csv(fh.name)
    fh.delete()

    fh = TempFile()
    m1.to_sif(fh.name)
    fh.delete()

    # m1 and m2 are identical. Adding them gets rid of duplicates so it should be
    # equal to itself.
    m1 == m1 + m2
    
