from cno.admin import install
import nose
import logging
import sys


def test_install():
    install.install_all_cellnopt_dependencies(verbose=True)



