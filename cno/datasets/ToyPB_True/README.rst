:Description: The PKN network **PKN-ToyPB_True.sif** is used in the reference below.
    It was generated manually based on the network contained in Saez-Rodriguez 
    et al. (2011). This network is designed to output a variety of dynamics (oscillations, 
    transience etc.) to test all CNOR formalisms.

Note that in the paper cited here below (MacNamara), the data was generated with
this model and CNORode software. A version of the model was used as a PKN. 
The differences are as follows:

========================================= =====================================================================
PKN-ToyPB                                 PKN-ToyPB_True
========================================= =====================================================================
map3k7-->mkk4 OR map3k7-->mkk4            This is a logical AND gate only
egfr-->sos OR `ph--|sos`                  This is a logical AND gate only (egfr-->sos AND ph--|sos )
connection from pi3k to rac added
connection from tnfr to pi3k added
new node (rac) between pi3k and map3k1
Node ask1 removed
========================================= =====================================================================

:Note: the MIDAS file MD-ToyPB_True.csv is identical to MD-ToyPB.csv

.. image:: https://github.com/cellnopt/cellnopt/blob/master/cno/datasets/ToyPB_True/ToyPB_True.png
   :alt: ToyPB_True figure



References
--------------

.. [1] **MacNamara et al.** 
    *State-time spectrum of signal transduction logic models* 
    Physical Biology, submitted (2012)

