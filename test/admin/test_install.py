from cno.admin import install
import nose
import logging
import sys

from nose.plugins.attrib import attr


@attr('skip_travis')
def test_install():
    install.install_all_cellnopt_dependencies(verbose=True)



