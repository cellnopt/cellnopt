"""This test file cover the cno.midas module"""
from nose.plugins.attrib import attr

from cno.io.midas import  XMIDAS, Trend
from cno.io.measurements import MIDASBuilder
from cno import cnodata

from easydev.easytest import assert_list_almost_equal, TempFile

from cno.testing import getdata
from cno.misc import CNOError

import numpy as np

filenames = getdata(pattern="MD*")


# Trying
def test_reading_and_saving():
    for filename in filenames:
        if filename.endswith('xml'):
            continue
        yield reading_and_saving, filename

def reading_and_saving(filename):
    # if multiple cell line, we must provide the cell line name
    if 'filtering' in filename:
        return

    if 'multiple' in filename and "wrong" not in filename:
        m1 = XMIDAS(filename, cellLine='C1')
        try:
            m1 = XMIDAS(filename)
            assert False
        except:
            assert True
    # if wrong MIDAS, a CNOError should be raised
    elif 'wrong' in filename:
        try:
            m1 = XMIDAS(filename)
            assert False
        except CNOError:
            assert True
            return
        except:
            assert False
    else:
        m1 = XMIDAS(filename)

    # we can also test the following: writing the XMIDAS into a file (MIDAS
    # format) and read it back. the 2 objects must be identical !
    f = TempFile()
    m1.to_midas(f.name)
    m2 = XMIDAS(f.name)
    f.delete()
    assert m1 == m2

    print(m1)

def test_read_xmidas_nofilename():
    m = XMIDAS()
    m.cellLine

def test_no_stimuli():
    m = XMIDAS(getdata("MD-test_one_stimulus.csv"))
    m.remove_stimuli('egf')
    m.names_stimuli

def test_no_inhibitors():
    m = XMIDAS(getdata("MD-test_one_stimulus.csv"))
    m.remove_inhibitors(['pi3k', 'raf1'])
    m.names_inhibitors

def test_incompatible_cellline():
    # multiple cell lines
    try:
        m = XMIDAS(getdata("MD-test_wrong_multiple_cellline.csv"), 'C1')
        assert False
    except CNOError:
        assert True

def test_xmidas():
    # constructor
    m = XMIDAS(cnodata("MD-ToyPB.csv"))
    m2 = XMIDAS(m)
    assert m == m2

    # accessors
    m.df
    assert len(m.names_cues)==4
    assert len(m.experiments) == 10
    assert len(m.names_stimuli)==2
    assert len(m.stimuli)==10
    assert len(m.inhibitors) == 10
    assert len(m.names_inhibitors) == 2
    assert m.cellLine == 'Cell'
    assert len(m.cellLines) == 1
    assert len(m.times ) == 16

    # getter
    assert m['p38'].shape == (160,)
    assert m['Cell','experiment_0', 0].shape == (6,)
    try:
        m.cellLine = "wrong"
        assert False
    except:
        assert True


    # ----------------------------------------- remove methods
    try:
        m.remove_species("dummy")
        assert False
    except:
        assert True
    m.remove_species("erk")
    assert len(m.species) == 5
    m.remove_times(18)
    m.remove_times([12,14])
    m.remove_experiments("experiment_0")
    m.remove_stimuli("p38")
    m.remove_inhibitors("pi3k")

    m.remove_inhibitors("dummy")
    m.remove_stimuli("dummy")

    mm = m.copy()
    try:
        mm.remove_cellLine("Cell")
        assert False
    except:
        assert True
    del mm

    # ------------------------------------------- rename methods
    m.rename_cellLine({"Cell": "HT29"})
    m.rename_species({"ap1":'AP1'})
    m.rename_inhibitors({'raf1':'RAF1'})
    m.rename_stimuli({'egf':'EGF'})
    m.rename_time({0:0,2:2*60,4:4*60})

    # others
    m.get_residual_errors()

def test_xmidas_average():
    m = XMIDAS(getdata("MD-test_average.csv"))
    m.average_replicates()
    m.average_replicates(inplace=True)

def test_operators():
    m = XMIDAS(cnodata("MD-ToyPB.csv"))
    # test *, +, and / operators
    assert (m + m) * 0.25 == m / 2
    assert (m - 0.5) + 0.5 == m
    assert (m - m).df.sum().sum() == 0


#plotting

def test_xmidas_corr():
    m = XMIDAS(cnodata("MD-ToyPB.csv"))
    m.corr()

def test_xmidas_plot():
    m = XMIDAS(cnodata("MD-ToyPB.csv"))
    m.create_random_simulation()
    m.plot(mode="trend")
    m.plot(mode="mse", vmax=.9) # vmax useful for mse only
    m.xplot(bbox=True)

    #m = XMIDAS(getdata("MD-unnorm.csv"))
    #m.plot(mode="mse", logx=True)
    #m.df -= 2
    ## TODO should raise an erro
    # m.plot(mode="mse")


def test_hcluster():
    m = XMIDAS(cnodata("MD-ToyPB.csv"))
    m.hcluster()

def test_boxplot():
    m = XMIDAS(cnodata("MD-ToyPB.csv"))
    m.boxplot(mode="time")
    m.boxplot(mode="experiment")

def test_radviz():
    m = XMIDAS(cnodata("MD-ToyPB.csv"))
    m.radviz()

def test_dicretise():
    m = XMIDAS(cnodata("MD-ToyPB.csv"))
    m.discretise(inplace=False)
    m.discretize(inplace=False)
    m.discretise(N=3)
    assert set(m.df.values.flatten()) == {0,1,.5}

def test_heatmap():
    m = XMIDAS(cnodata("MD-ToyPB.csv"))
    m.heatmap()



# scaling
def test_scaling():
    m = XMIDAS(cnodata("MD-ToyMMB.csv"))
    #m.scale_max_by_experiments()
    #assert (m.df == 1).sum().sum() == 40

    m = XMIDAS(m)
    m.scale_max()
    #assert (m.df == 1).sum().sum() == 40

    # SCALE MIN/MAX
    m1 = XMIDAS(cnodata("MD-ToyMMB.csv"))
    m2 = XMIDAS(cnodata("MD-ToyMMB.csv"))
    m2.scale_min_max()
    # because min is 0 and max is 1, scale_min_max should not affect the datafame
    assert m1 == m2


    #m = XMIDAS(m)
    #m.scale_max_by_experiments()

    #m = XMIDAS(m)
    #m.scale_min_max_by_experiments()

    m = XMIDAS(m)
    m.scale_max_across_experiments()

    m = XMIDAS(m)
    m.scale_min_max_across_experiments()

    # same but with inplace set to False
    m = XMIDAS(m)
    m.scale_max(inplace=False)
    m = XMIDAS(m)
    m.scale_min_max(inplace=False)
    #m = XMIDAS(m)
    #m.scale_max_by_experiments(inplace=False)
    #m = XMIDAS(m)
    #m.scale_min_max_by_experiments(inplace=False)
    m = XMIDAS(m)
    m.scale_max_across_experiments(inplace=False)
    m = XMIDAS(m)
    m.scale_min_max_across_experiments(inplace=False)


def test_shuffle():
    m = XMIDAS(cnodata("MD-ToyPB.csv"))

    # shuffling over signals should keep the sum over signal constan
    a1 = m.df.sum(axis=0)
    for this in ['signal', 'column']:
        m.shuffle(mode=this)
        a2 = m.df.sum(axis=0)
        assert_list_almost_equal(a1, a2)
        m.reset()

    # shuffling over index keep sum over a given experiment constant
    a1 = m.df.sum(level="experiment").sum(axis=1)
    m.shuffle(mode="index")
    a2 = m.df.sum(level="experiment").sum(axis=1)
    assert_list_almost_equal(a1, a2)
    m.reset()

    m.shuffle(mode="timeseries")

    m.reset()
    a1 = m.df.sum().sum()
    m.shuffle(mode="all")
    a2 = m.df.sum().sum()
    assert_list_almost_equal([a1], [a2])


# NORMALISATION


def test_norm_time():
    m = XMIDAS(getdata("MD-test_unnorm.csv"))
    m.normalise(mode="time")

def test_norm_control():
    m = XMIDAS(getdata("MD-test_unnorm_exp.csv"))
    m.normalise(mode="control")


#MEASUREMENTS

def test_xmidas_experiments():
    m = XMIDAS(cnodata("MD-ToyPCB.csv"))
    ms = m.to_measurements()
    assert len(ms) == 60
    # we can create a MIDASBuilder, exoprts and can check we get the same ?
    # does not work for now
    mb = MIDASBuilder()
    mb.add_measurements(ms)
    m2 = mb.xmidas
    assert m == m2


# TREND

def test_trend():
    xm = XMIDAS(cnodata("MD-ToyPB.csv"))
    ts = xm.df['ap1']['Cell']['experiment_0']
    # set a trend instance
    trend = Trend()
    trend.set(ts)
    trend.plot()
    trend.get_bestfit_color()



def test_constructor_time():
    m = XMIDAS(getdata("MD-test_time0_to_duplicate.csv"))
    assert m.times == [0,30,180]
    assert m.df.shape == (12,1)
    assert np.all(m.experiments.values == np.array([[ 0.,  0.,  0.,  0.],
       [ 1.,  0.,  0.,  0.],
       [ 0.,  0.,  1.,  0.],
       [ 1.,  0.,  1.,  0.]]))

    assert np.all(m.df.values == np.array([[350],
       [409],
       [290],
       [318],
       [377],
       [258],
       [350],
       [560],
       [344],
       [318],
       [485],
       [269]]))

def test_constructor_no_time_zero():
    m = XMIDAS(getdata("MD-test_no_time0.csv"))

def test_filtering():
    m = XMIDAS(getdata("MD-test_filtering.csv"), exclude_rows={'ID:type':'blank'})
    assert m.df.sum().sum() == 2072


"""


def _test_compare_pycno_vs_cnor():
    # test all sif files
    for f in filenames:
        yield compare_pycno_vs_cnor, f


def _compare_pycno_vs_cnor(filename):
    #check that the sif file read by CNOR and pyCNO are identical
    filename = cnodata(filename)
    try:
        from cellnopt.wrapper import rpack_CNOR as cnor
        from cellnopt.core import midas

        mr = cnor.readMIDAS(filename)
        m = midas.MIDASReader(filename)

        assert sorted(list(mr.rx2('DAcol'))) == sorted([x+1 for x in list(m.DAcol)])
        assert sorted(list(mr.rx2('DVcol'))) == sorted([x+1 for x in list(m.DVcol)])
        assert sorted(list(mr.rx2('TRcol'))) == sorted([x+1 for x in list(m.TRcol)])
    except:
        pass

    # TODO(tc): seems that we need to transpose one matrix. do not seem neat
    #assert numpy.array(m.dataMatrix).all() == numpy.array(mr.dataMatrix).transpose().all()

"""
