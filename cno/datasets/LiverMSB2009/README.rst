:Description:  
:Model: **PKN-LiverMSB2009.sif** provides a PKN used in the reference here below and **PKN-LiverMSB2009.xml** is the
    SBML-qual version.
:Data: data in MIDAS format is in **MD-LiverMSB2009.csv**


.. image:: https://github.com/cellnopt/cellnopt/blob/master/cno/datasets/LiverMSB2009/LiverMSB2009.png
   :alt: LiverMSB2009 figure
   :scale: 30%

-

.. image:: https://cdn.rawgit.com/cellnopt/cellnopt/master/cno/datasets/LiverMSB2009/LiverMSB2009.svg
   :alt: LiverMSB2009 figure
   :target: 
   
-
    
references
----------------


.. [1] **J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger.**
   *Discrete logic modeling as a means to link protein signaling networks with functional analysis of mammalian signal transduction*
   Molecular Systems Biology, 5:331, 2009

Notes
--------

There are differences with the version provided in the journal because PKN and MIDAS were not compatible. Here are the differences. In the MIDAS file:
    
    - add :CellLine after TR:HepG2 
    - EGFRi-Gefin --> EGFRi
    
    In the PKN::

        erk12 to ERK
        MEK12 to MEK
        JNK12 to JNK 
        IRS1S to IRS-1S

    * Same modifications for the PKN. In addition, the following reactions are
      removed from the PKN provided on the journal website::

        STAT6   1   dummy
        CASP37  1   dummy
        BCA 1   dummy







