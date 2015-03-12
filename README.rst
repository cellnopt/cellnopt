cellnopt/cno
=============

.. image:: https://badge.fury.io/py/cellnopt.svg
    :target: https://pypi.python.org/pypi/cno

.. image:: https://pypip.in/d/cellnopt/badge.png
    :target: https://crate.io/packages/cellnopt

.. image:: https://secure.travis-ci.org/cellnopt/cellnopt.png
    :target: http://travis-ci.org/cellnopt/cellnopt

.. image:: https://coveralls.io/repos/cellnopt/cellnopt/badge.png?branch=master 
   :target: https://coveralls.io/r/cellnopt/cellnopt?branch=master 

.. image:: https://landscape.io/github/cellnopt/cellnopt/master/landscape.png
   :target: https://landscape.io/github/cellnopt/cellnopt/master

.. image:: https://badge.waffle.io/cellnopt/cellnopt.png?label=ready&title=Ready 
   :target: https://waffle.io/cellnopt/cellnopt

:Note: Cellnopt is tested under Python2.7 and portage to Python3 is on the way
:Contributions: Please join `cellnopt github project <https://github.com/cellnopt/cellnopt>`_ and share your
                 `notebooks <https://github.com/cellnopt/cellnopt/tree/master/notebooks>`_
:Issues: Please use `cellnopt issues and bug report <https://github.com/cellnopt/cellnopt/issues>`_
:Documentation: Some documentation (mostly API) are also available on `pypi <http://pythonhosted.org//cno/>`_

This is a stable version that contains datasets and tools to manipulate
models and data sets used in CellNOpt project. The sub-package cno.io is ready and 
is a replacement for `cellnopt.core <https://pypi.python.org/pypi/cellnopt.core>`_ package.
cno.datasets sub-package is a replacement for `cellnopt.data <https://pypi.python.org/pypi/cellnopt.data>`_ package.

As for the formalisms, pipelines are provided for each formalism and are available in the 
boolean, ode and fuzzy packages. 

:See also: http://www.cellnopt.org for the context.

Installation
===============

From Pypi / pip::

    sudo pip install cno
    
From github::

    # Get the source either from http
    git clone http://github.com/cellnopt/cellnopt.git
    # or SSH
    # git clone git@github.com:cellnopt/cellnopt.git
    cd cellnopt
    sudo python setup.py install






