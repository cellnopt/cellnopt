from cno.io import ADJ2SIF, SIF
from cno.testing import getdata
# FIXMED import CNOGraph
import os
from  tempfile import mkstemp

def test_adj2sif():

    f1 = getdata("test_adjacency_matrix.csv")
    f2 = getdata("test_adjacency_names.csv")

    s = ADJ2SIF(f1, f2)
    s.to_sif()

    fd, name = mkstemp(suffix=".sif")
    s.to_sif(name)
    s2 = SIF(name)
    s2.nodes1 == ["A", "A"]
    s2.nodes2 == ["B", "C"]
    s2.edges == ["1", "1"]

    # FIXME
    #c = CNOGraph(s.G)

def test_load():
    a = ADJ2SIF()
    f1 = getdata("test_adjacency_matrix.csv")
    f2 = getdata("test_adjacency_names.csv")
    a.load_adjacency(f1)
    a.load_names(f2)

    # FIXME
    # CNOGraph(a.G) # just to check that it can be imported via cnograph

