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

There is one restriction: **do not use dash in the name** (e.g., Toy-MMB) because automatic import in Python may fail.

Otherwise please follow those rules:

  #. `upper camel case convention <http://en.wikipedia.org/wiki/CamelCase>`_.
  #. A directory that contain a variant of an existing directory should
     reuse the identifier of the existing directory. Then, append a second identifier
     separated by an underscore or no separation.

:Examples: ::

       ToyExample
       
       ToyExample_Feedback or ToyExampleFeedback
       
       ToyExample_Feedback_bis

The **ToyExample_Feedback_bis** is therefore based on **ToyExample_Feedback** that is itself
a variant of **ToyExample**. In this situation, we could have called the third
directory simply **ToyExample_FeedbackBis** since it is also a variant of **ToyExample**.


2. Contents
---------------

In a directory, you must provide 3 files:

    #. a network (SIF format) named **PKN-Identifier.sif**
    #. a README.rst file (contents does not require any specific template)
    #. an __init__.py file
    
Optional data files also be provided:

    #. a MIDAS file (csv format) named **MD-Identifier.csv**
    #. a network in SBML qual format (extension must be .xml)


1.1 PKN file
~~~~~~~~~~~~~~
The PKN file must be in SIF format and named **PKN-<identifier>.sif**

1.2 MIDAS file
~~~~~~~~~~~~~~~
The data file must be MIDAS format and named **MD-<identifier>.csv**
If not provided because it already exist in another directory, use the __init__ file (see below)

1.3 README file
~~~~~~~~~~~~~~~~~~~~
Must be in RST format (like this document). Please, see examples within this github repository.

1.4 __init__ file
~~~~~~~~~~~~~~~~~~~~~~

The __init__ file can be empty, in which case the name of the SIF, MIDAS and SBSML files will be inferred
from the directory name. If not MIDAS is provided or for some reasons the SIF (and SBML) are named differently, 
you can use a metadata dictionary into the __init__ file:

```

metadata = {
    'data': '../ToyMMB/MD-ToyMMB.csv'
}    

```


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







