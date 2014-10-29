

#############################
CellNOpt documentation
#############################

.. topic:: cellnopt

    .. image:: network.png
        :width: 30%


Motivation 
###########


**Cellnopt** is a pure Python library that provides core functions to manipulate signalling data. It provides functions to create and manipulate networks that can be used later on with MIDAS files (phosphorylation measurements). You can read/write data formats such as :ref:`sif` and :ref:`midas` that are used in `CellNOptR <http://bioconductor.org/packages/release/bioc/html/CellNOptR.html>`_.

It does not provide any logical formalism to simulate or optimise the network to
the data. It is rather intended to be a package dedicated to the pre and post
processing of signalling pathways. 


Installation
###################

**cellnopt.core** has a dependencies on **pandas**, which requires cython. So, you need to install cython first. 

All other dependencies should be installed automatically (networkx, numpy, matplotlib).

::

    pip install easydev
    pip install cython
    pip install cellnopt



User guide
##################


.. toctree::
    :maxdepth: 2
    :numbered:

    quickstart.rst
    userguide.rst
    applications.rst

References
##################


.. toctree::
    :maxdepth: 2
    :numbered:

    references


.. toctree::

    biblio.rst


.. toctree::

    todo.rst
    ChangeLog.rst
    authors.rst
    contribute
    faqs.rst
