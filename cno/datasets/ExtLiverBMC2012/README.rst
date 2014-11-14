:Description:  **ExtLiverBMC2012** derives from :ref:`ExtLiverPCB`. The data is the same so no
new data file is provided. The PKN model is slightly different (See below).
:Model: The model to be used **PKN-ExtLiverBMC2012.sif**. Note, However, that the directory 
also contains **PKN-ExtLiverMSBmodUP.sif**, which has the same topology but names are UniProt identifiers. 
:Data:
    **MD-ExtLiverBMC2012.csv** is exactly the same file as **MD-ExtLiverPCB.csv**. It is normalised.

.. warning:: Note that they are a few variants of this data/PKN.


.. image:: https://github.com/cellnopt/cellnopt/blob/master/cno/datasets/ExtLiverBMC2012/ExtLiverBMC2012.png

References
----------------

.. [1] **Terfve, C. and Cokelaer, T. and Henriques D. and MacNamara A. and Gonçalves E. and Morris K.M. and van Iersel M. and Lauffenburger A.D. and Saez-Rodriguez J.**,
    *CellNOptR: a flexible toolkit to train protein signaling networks to data using multiple logic formalisms*,
    BMC Systems Biology 2012, 6:133,
    `Citation <http://www.biomedcentral.com/1752-0509/6/133/abstract>`_


The second set of files are the original from the reference and differ from the
first set only by their naming convention being uniprot.

Notes
-------------
Differences between **PKN-ExtLiverBMS2012.sif** with :ref:`ExtLiverPCB` model::

    * 3 additional links::

        erk12   -1  sos
        AKT     -1  ras
        IRS1s   -1  irs1t

    * One feedback link removed::

        IL6r    1   PI3K

Other files available
--------------------------

    * MCP_HepG2mod4.csv (UNNORMALISED!!)
    * MCP_PriHumod5.csv (UNNORMALISED!!). Same as in :ref:`ExtLiverPriHu-MCP2010`
