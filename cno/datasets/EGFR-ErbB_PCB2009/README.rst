:Description: **EGFR-ErbB_PCB2009**
:Model: The network was generated from the *reactions* file provided by the author and parsed by **Reactions** (see T.C.). 
:Data: the data is available in **MD-EGFR-ErbB_PCB2009.csv** and is a copy of *./data/BoolData-Exp1.csv*. Other MIDAS files can be found in ./data and were provided by the authors. Note that the header was edited to fix the cell line  tag.


.. image:: https://github.com/cellnopt/cellnopt/blob/master/cno/datasets/EGFR-ErbB_PCB2009/EGFR-ErbB_PCB2009.png

```
del c.graph_options['graph']['ratio']
c.graph_options['graph']['ranksep'] = 0.2
c.graph_options['graph']['size']= '15,15'
c.plot()
```

References
------------

.. [1] **R. Samaga, J. Saez-Rodriguez, L. Alexopoulos, P. K. Sorger, S. Klamt.**, 
   *The logic of EGR/ErbB signaling: theoretical properties and analysis of high-throughput data.* 
   PLoS Comp. Biol., 5(8): e1000438, 2009.
   `Citation <http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000438>`_

Notes
-------
The network file has also been modified to be compatible with the MIDAS file.
Species that are measured have the same name as the ones found in the MIDAS file
taking care of lower and upper cases. 

NOLIG changed into CYTO.
