**EGFR-ErbB_PCB2009**

This model is used in the references here below. The network was generated from the *reactions* file provided by the author and parsed by **read_reactions.Reactions** (see T.C.). The data (MIDAS/csv file) is called *EGFR-ErbB_PCB2009.csv* and is a copy of *./data/BoolData-Exp1.csv*. Other MIDAS files can be found in ./data and were provided by the authors. Note that the header of *EGFR-ErbB_PCB2009.csv* is different from the original file::

    TR:HepG2:Cell type 

is now::

    TR:HepG2:CellLine 

The network file has also been modified to be compatible with the MIDAS file.
Species that are measured have the same name as the ones found in the MIDAS file
taking of of lower and upper cases. 

NOLIG changed into CYTO.

:References:

.. [1] **R. Samaga, J. Saez-Rodriguez, L. Alexopoulos, P. K. Sorger, S. Klamt.**, 
   *The logic of EGR/ErbB signaling: theoretical properties and analysis of high-throughput data.* 
   PLoS Comp. Biol., 5(8): e1000438, 2009.
   `Citation <http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000438>`_
