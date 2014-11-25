from cno.admin import DistributeRPackage
from nose.plugins.attrib import attr


@attr('skip_travis')
def test_distribute():

    dp = DistributeRPackage('CellNOptR')
    dp.distribute()
    dp.clean()

    try:
        dp = DistributeRPackage('MEIGOR')
        dp.distribute()
        assert False
    except:
        assert True


    dp.help()


# TODO
def test_distribute_script():
    pass


# FIXME WONT DO
def install_all_cellnopt_dep():
    pass
