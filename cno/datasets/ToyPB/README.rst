:Description: The PKN network **PKN-ToyPB.sif** is used in the reference below.
    It was generated manually based on the network contained in Saez-Rodriguez 
    et al. (2011). This network is designed to output a variety of dynamics (oscillations, 
    transience etc.) to test all CNOR formalisms.

Note that in the paper cited here below (MacNamara), the data was generated
using CNORode software with a modified version of the model. The **True** Model
is called **PKN-ToyPB_True.sif**. The differences are as follows:

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


.. image:: https://github.com/cellnopt/cellnopt/blob/master/cno/datasets/ToyPB/PKN-ToyPB.png
   :width: 200pt
   :height: 100pt
   :align: center
   :alt: ToyPB figure


- SBML format available: cnodata("PKN-ToyPB.xml")


References
--------------

.. [1] **MacNamara et al.** 
    *State-time spectrum of signal transduction logic models* 
    Physical Biology, submitted (2012)

.. note:: The directory contains a SIF file named ToyModelPB.sif but also a
    PKN-ToyPB_shuffled.sif that has exactly the same topology. However, this second
    file has a different arrangements of reactions. When compressed and expanded, 
    the resulting  models are different. This is a known bug and these 2 files can 
    be used to debug that issue. The file MD-ToyPB_2timePoints.csv is a subset of 
    ToyModelPB.csv that contains the time points 0, 10 and 30 seconds only.
