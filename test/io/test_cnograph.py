from cno.io import CNOGraph, XMIDAS, SIF
from cno import cnodata, getdata
import tempfile


sif = getdata("PKN-Test.sif")
midas = getdata("MD-Test.csv")


def test_class_constructor():
    sif = SIF(cnodata("PKN-ToyPB.sif"))
    midas = XMIDAS(cnodata("MD-ToyPB.csv"))
    c = CNOGraph(sif, midas)
    N = len(c)
    c.compress()
    assert len(c) == 12
    c._decompress()
    # assert len(c) == N
    c.plotdot2()
  

def test_plotting():
    c = CNOGraph(sif, midas)
    c.plotdot(legend=True, show=False)
    c.centrality_closeness()
    c.plotdot(node_attribute="centrality_closeness", show=False) 
 

def test_plotFeedback():
    c = CNOGraph(sif, midas)
    c.plotFeedbackLoopsSpecies()
    c.plotFeedbackLoopsHistogram()

def test_summary():
    c = CNOGraph(sif, midas)
    #c.summary()

def test_class2():
    G = CNOGraph()
    G.add_node("a")
    G.add_node("b")
    G.add_node("c")
    G.add_node("d")
    G.add_edge("a", "b", link="+")
    G.remove_node("c")
    G["a"]


def test_add_cycle():
    c = CNOGraph()
    c.add_cycle(["A", "B", "C"], link="+")
    assert len(c.nodes()) == 3
    assert len(c.edges()) == 3

def test_operators():
    c1 = CNOGraph()
    c1.add_edge("A","B",link="+")
    c1.add_edge("A","C",link="+")
    c2 = CNOGraph()
    c2.add_edge("A","E", link="+")
    c2.add_edge("C","E", link="+")
    c2.plotdot(show=False)
    assert c1 != c2

    c_add = c1+c2
    assert c_add == c1+c2

    c_radd = c1
    c_radd += c2

    assert c_add == c_radd

    c3 = CNOGraph(); 
    c3.add_edge("A", "B", link="+")
    c_sub = c_add -c3
    sorted(c_sub.nodes()) == ["C", "E"]


    c4 = CNOGraph()
    c4.add_edge("A+B" ,"C", link="+")

def test_clean_orphan_ands():
    c = CNOGraph()
    c.add_edge("A", "B", link="+")
    c.add_edge("A", "C", link="+")
    c.add_edge("Aprime", "C", link="+")
    c.add_edge("C", "D", link="+")
    c.add_edge("B", "D", link="+")
    c.expand_and_gates(2)
    c.remove_node("C")
    c.remove_node("Aprime")
    assert len(c) == 3


def test_expandORgates():
    c = CNOGraph()
    c.add_edge("A", "C", link="-")
    c.add_edge("B", "C", link="+")
    c.expand_and_gates()
    c.remove_edge("B", "C")

    assert len(c.edges()) == 4
    c.expand_or_gates()
    assert len(c.edges()) == 5



def test_centrality():
    c = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
    c.centrality_betweeness()
    c.centrality_degree()
    c.degree_histogram(show=True)
    c.draw()

def test_preprocessing():
    # functional test of the preprocessing steps
    c = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
    c.cutnonc()
    c.compress()
    c.expand_and_gates()

    c = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
    c.preprocessing()

    c = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
    c.cutnonc()
    c.compress()
    c.expand_and_gates(2)
    assert len(c) == 17
    c.expand_and_gates(3) # there is no nodes with 3 inputs so length is the same
    assert len(c) == 17 

    # remove and nodes could be call several times without effect.
    c.remove_and_gates()
    assert len(c) == 12


def test_export():
    fh = tempfile.NamedTemporaryFile()
    c = CNOGraph(sif, midas)
    c.export2sif(fh.name)
    c.export2SBMLQual(fh.name)
    c.export2gexf(fh.name)

def test_others():
    import os
    c = CNOGraph(sif, midas)
    c.adjacencyMatrix()
    c.export2json("test.json")
    c.loadjson("test.json")
    os.remove("test.json")
    print(c)
    c.reacID
    c.namesSpecies
    c.dependencyMatrix()
    c.plotAdjacencyMatrix()
    c.get_same_rank()
    c.dot_mode = "signals_bottom"
    c.get_same_rank()
    c.lookfor("akt")
    c.lookfor("EGFR")


def test_set_operators():

    c1 = CNOGraph()
    c1.add_edge("A","C",link="+")
    c1.add_edge("A","B",link="+")
    c2 = CNOGraph()
    c2.add_edge("A", "C", link="+")
    c1.intersect(c2).plotdot(show=False)
    c3 = c1.intersect(c2)
    assert sorted(c3.nodes()) == ['A', 'C']


    c1 = CNOGraph()
    c1.add_edge("A","C",link="+")
    c1.add_edge("A","B",link="+")
    c2 = CNOGraph()
    c2.add_edge("A", "C", link="+")
    c2.add_edge("A", "D", link="+")
    c3 = c1.union(c2)
    assert sorted(c3.nodes()) == ['A', 'B', 'C', 'D']
    assert len(c3.edges()) == 3

    c3 = c1.difference(c2)
    assert c3.nodes() == ['B']




def test_merge_split():
    c = CNOGraph()
    c.add_edge("AKT2", "B", link="+")
    c.add_edge("AKT1", "B", link="+")
    c.add_edge("A", "AKT2", link="+")
    c.add_edge("A", "AKT1", link="+")
    c.add_edge("C", "AKT1", link="+")
    c.merge_nodes(["AKT1", "AKT2"], "AKT")
    assert sorted(c.nodes()) == ["A", "AKT", "B", "C"]

    c = CNOGraph()
    c.add_edge("AKT2", "B", link="+")
    c.add_edge("AKT1", "B", link="+")
    c.add_edge("A", "AKT2", link="+")
    c.add_edge("A", "AKT1", link="-")
    c.add_edge("C", "AKT1", link="+")
    c.merge_nodes(["AKT1", "AKT2"], "AKT")

    assert sorted(c.nodes()) == ["A", "AKT", "B", "C"]

    c.split_node("AKT", ["AKT1", "AKT2"])


def test_sif():
    c = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyPB.csv"))
    c2 = CNOGraph(c.sif())
    assert c == c2

def test_import_sif_with_and_gates():
    fh = tempfile.NamedTemporaryFile(delete=False)
    s = SIF()
    s.add_reaction("A^B=C")
    c = CNOGraph(s)
    c.export2sif(fh.name)
    c2 = CNOGraph(fh.name)
    assert (c == c2) == True
    fh.delete = True
    fh.close()

def test_1cue_noinh():
    s = get_share_file("PKN-test_1cue_noinhibitor.sif")
    m = get_share_file("MD-test_1cue_noinhibitor.csv")
    c = CNOGraph(s, m)


def test_expand_4_and_gates():
    s = get_share_file("PKN-test_4andgates.sif")
    m = get_share_file("MD-test_4andgates.csv")

    c2 = CNOGraph(s,m)
    c1 = CNOGraph(s,m)
    c1.expand_and_gates()
    c1.remove_and_gates()

    assert c1 == c2

def test_diamond():
    s = get_share_file("PKN-test_diamond.sif")
    m = get_share_file("MD-test_diamond.csv")
    c = CNOGraph(s,m)
    assert len(c.nodes()) == 4

    c.preprocessing()
    #TODO: check that there is one node not compressed due to the 2 links 1 stimuli one inhibitor

    
