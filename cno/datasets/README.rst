Description
###############

The **cno.datasets** package gathers published pairs of model/data files. You can add new directories as explained here below. You can browse the `source directory <https://github.com/cellnopt/cellnopt/tree/master/cno/datasets>`_ to visualise the models. Each sub-directory is actually a Python package, which means you can easily retrieve information, pathnames and visualise the model. For instance the **ToyMMB** can be retrieved as follows

.. doctest::

    >>> from cno.datasets import ToyMMB
    >>> ToyMMB.model_filename
    cno/datasets/ToyMMB/PKN-ToyMMB.sif
    >>> cno.datasets.ToyMMB.name
    ToyMMB
    >>> cno.datasets.ToyMMB.plot()


How to add a new directory/package ?
=====================================

1. Decide on an identifier (and naming convention)
----------------------------------------------------

If you want to add a model/data, first figure out a unique identifier, which will be used to 
create a new sub-directory.

There is no strict rules on naming the identifier but please use 

  #. `upper camel case convention <http://en.wikipedia.org/wiki/CamelCase>`_.
  #. A directory that contain some data based on an existing directory should
     reuse the identifier of the existing directory and add a second identifier
     separated by an underscore.

:Examples: ::

       ToyExample
       
       ToyExample_Feedback
       
       ToyExample_Feedback_bis

The **ToyExample_Feedback_bis** is therefore based on **ToyExample_Feedback** that is itself
a variant of **ToyExample**. In this situation, we could have called the third
directory simply **ToyExample_FeedbackBis** since it is also a variant of **ToyExample**.


2. Contents
---------------

In a directory, you must provide 4 files:

    #. a network (SIF format) named **PKN-Identifier.sif**
    #. a MIDAS file (csv format) named **MD-Identifier.csv**
    #. a README.rst file (contents does not require any specific template)
    #. an __init__.py file

1.1 PKN file
~~~~~~~~~~~~~~
The PKN file must be in SIF format and named **PKN-<identifier>.sif**

1.2 MIDAS file
~~~~~~~~~~~~~~~
The data file must be MIDAS format and named **MD-<identifier>.csv**

1.3 README file
~~~~~~~~~~~~~~~~~~~~
Must be in RST format (like this document). Please, see examples within this github repository.

1.4 __init__ file
~~~~~~~~~~~~~~~~~~~~~~

Just copy and paste this file into your directory changinf the name (i.e. here below, 
change ToyPB to your identifier)::

    import os
    from . import __path__

    name = 'ToyPB'

    __all__ = ['description', 'model_filename', 'data_filename', 'name']

    model_filename = __path__[0] + os.sep + "PKN-{0}.sif".format(name)
    data_filename = __path__[0] + os.sep + "MD-{0}.csv".format(name)
    description = open(__path__[0] + os.sep + "README.rst").read()

    __doc__ += description

    def plot():
        from cellnopt.core import CNOGraph
            CNOGraph(model_filename, data_filename).plot()

1.5 Others
~~~~~~~~~~~~~~~
You can of course add as much files as you want but please keep size as low as possible. 
Examples of other files that can be provided are: 

* SVG file of the network
* dot file of the network
* SBML-qual format

:Note: When using **cnodata** function provided in CellNOpt, all sub-directories with an 
__init__.py file are scanned automatically and the SIF, MIDAS and SBML-qual files with their full path 
are known. Type **cnodata()** without argument to get the list of
filenames that are available. Type **cnodata(identifier)** to retrieve the full path of the file (see below)

3. Usage
--------------
Once you have added a sub-directory, users and developers can then access to your data easily::

    from cno import cnodata
    pkn = cnodata("PKN-identifier.sif")
    midas = cnodata("MD-identifier.csv")







