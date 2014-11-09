Description
###############

The **cno.datasets** package gathers published pairs of model/data files. You can add new directories as explained here below. You can browse the `source directory <https://github.com/cellnopt/cellnopt/tree/master/cno/datasets>`_ to visualise the models. Each sub-directory is actually a Python package, which means you can easily retrieve information, pathnames and visualise the model. For instance the **ToyMMB** can be retrieved as follows:: 

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


All model and data sets are stored in the **share/data** directory within a
dedicated sub-directory. This sub-directory may contain several files but no
sub-directories. 


Although the directory names are not important because **cellnopt.data** scans the entire **share/data**
directory, nevertheless it should reflect its contents and will be used as an identifier for the filenames.
The naming conventions for the directories are simple:

  #. The directory should use  `upper camel case convention <http://en.wikipedia.org/wiki/CamelCase>`_.
  #. A directory that contain some data based on an existing directory should
     reuse the identifier of the existing directory and add a second identifier
     separated by an underscore.

.. topic:: Examples: 

    ::

       ToyExample
       ToyExample_Feedback
       ToyExample_Feedback_bis

The **ToyExample_Feedback_bis** is therefore based on **ToyExample_Feedback** that is itself
a variant of **ToyExample**. In this situation, we could have called the third
directory simply **ToyExample_FeedbackBis** since it is also a variant of **ToyExample**.

The choice is let to the user but note that the generic name is then used to name files contained in it. So, the option **ToyExample_FeedbackBis** may look nicer.



Contents
-----------

In a directory, you should provide at least 3 files:

    #. a network (SIF format) named **PKN-Identifier.sif**
    #. a MIDAS file (csv format) named **MD-Identifier.csv**
    #. a description file (see below) called **description.txt**

.. warning:: note that the model Prior Knowledge Network (PKN) must starts with **PKN-** and end with **.sif** and the MIDAS filenames must start with **MD-** and end with **.csv**.

You can of course add as much files as you want. Examples are: 

 * SVG file of the network
 * dot file of the network
 * data file using reactions/metabolites format like in CNA
 * file called Type.NA if your model contains AND gates. This file is an
   attribute file for Cytoscape and is also used in the matlab version of CNO to
   import AND gates


aliases
-------------

We propose the above convention to start PKN model with **PKN-** and MIDAS file with **MD-**. The reason being that 
quite a few models and data filenames were using various conventions (e.g., MIDAS file where tagged **Data** or **data** or **MIDAS** and could appear anywhere in the filenames). 

Yet, some people used to use some specific filenames. If so, you can complete a
special file called aliases.py that may be used in some specific tools such as
:meth:`cellnopt.data.cnodata` to fetch a model or data set.

The file looks like::

    aliases = {
        "PKN-ToyMMB.sif": ["ToyModelMMB.sif", "ToyModelMKM.sif"],
        "MD-ToyMMB.csv": ["ToyModelMMB.csv", "ToyModelMKM.csv"]}



The description file
------------------------

The description file is compulsary. It is used to build this documentation
automatically and is of course important to keep track of the papers it was
first published in, authorship and so on.

The layout of the description file should be as faithful (keeping space and
return carriage!!) as possible to the following layout so that it can be 
intrepreted by the documentation builder automatically::

    Description of the model, its origins and main features. Can be as long as
    you want.

    :References:

    **authors of a paper**
    *title of the paper*
    Reference
    `Citation <http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=XXXXXXXXX>`_

    **authors of a paper**
    *title of the paper*
    Reference
    `Citation <http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=XXXXXXXXX>`_

    .. note:: alias to PKN-ToyMMB.sif are: ToyModelMMB.sif,  ToyModelMKM.sif (optional)
    .. note:: another note (optional)



Filename convention
------------------------------------------------

MIDAS files should start with "MD-" followed by a tag and the extension must be
".csv"::


    MD-TAG1.csv

SIF files should start with "PKN-" followed by a tag and the extension must be
".sif"::

    PKN-TAG1.sif

TAG1 is a label corresponding to the model. Variant of a file should have second tag as follows::

    PKN-TAG1_TAG2.sif

A compressed and expanded model to be saved could be saved as follows::

    CEN-Tag1.sif




