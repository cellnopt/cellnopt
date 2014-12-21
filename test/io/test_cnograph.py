from cno.io import CNOGraph, XMIDAS, SIF
from cno import cnodata, getdata, CNOError
import tempfile
import random
import os


sif = getdata("PKN-test_cnograph.sif")
midas = getdata("MD-test.csv")


def test_class_constructor():
    sif = SIF(cnodata("PKN-ToyPB.sif"))
    midas = XMIDAS(cnodata("MD-ToyPB.csv"))
    c = CNOGraph(sif, midas)
    N = len(c)
    c.compress()
    assert len(c) == 12
 
    # check constructor with another cnograph
    c2 = CNOGraph(c)
    assert c == c2
    # check that this is a copy not a reference
    c.remove_node('egf')
    assert c2 != c

    # check sbmlqual constructor
    c = CNOGraph(getdata("PKN-test_sbmlqual.xml"))
    assert "A^B=C" in c.reactions
    assert len(c.nodes()) == 5

    # test empty one
    CNOGraph()
    


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
    c2.plot(show=False, cmap='copper')
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

    c5 = CNOGraph()
    c5.add_edge("A+!B", "C", link="+")
    assert c5.reactions == sorted(['A=C', '!B=C'])



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


def test_to_sif():
    for filename in ['ToyPB', 'LiverDREAM', 'ToyMMB']:
        yield to_sif, filename


def to_sif(filename):
    c = CNOGraph(cnodata("PKN-"+filename+'.sif'))
    assert c == CNOGraph(c.to_sif())


def test_set_operators():
    c1 = CNOGraph()
    c1.add_edge("A","C",link="+")
    c1.add_edge("A","B",link="-")
    c2 = CNOGraph()
    c2.add_edge("A", "C", link="+")
    c1.intersect(c2).plot(show=False, cmap='green')
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


def test_import_sif_with_and_gates():
    fh = tempfile.NamedTemporaryFile(delete=False, suffix=".sif")
    s = SIF()
    s.add_reaction("A^B=C")
    c = CNOGraph(s)
    c.to_sif(fh.name)
    print(fh.name)
    c2 = CNOGraph(fh.name)
    print(c)
    print(c2)
    assert c == c2 
    fh.delete = True
    fh.close()

    # test incorrect extension
    fh = tempfile.NamedTemporaryFile(delete=False)
    c.to_sif(fh.name)
    try:
        CNOGraph(fh.name)
        assert Fales
    except:
        assert True
    fh.delete = True
    fh.close()


def test_1cue_noinh():
    s = getdata("PKN-test_1cue_noinhibitor.sif")
    m = getdata("MD-test_1cue_noinhibitor.csv")
    c = CNOGraph(s, m)


def test_expand_4_and_gates():
    s = getdata("PKN-test_4andgates.sif")
    m = getdata("MD-test_4andgates.csv")

    c2 = CNOGraph(s,m)
    c1 = CNOGraph(s,m)
    c1.expand_and_gates()
    c1.remove_and_gates()
    assert c1 == c2

def test_diamond():
    s = getdata("PKN-test_diamond.sif")
    m = getdata("MD-test_diamond.csv")
    c = CNOGraph(s,m)
    assert len(c.nodes()) == 4
    c.preprocessing()
    #TODO: check that there is one node not compressed due to the 2 links 1 stimuli one inhibitor


def test_check_compatible_midas():
    #try:
    #    c = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyMMB.csv"))
    #    assert False
    #except CNOError:
    #    assert True
    c = CNOGraph(cnodata("PKN-ToyPB.sif"), cnodata("MD-ToyMMB.csv"))

def test_renaming_of_ang_gates_after_compression():
    c = CNOGraph(cnodata("PKN-ToyPB_True.sif"), cnodata("MD-ToyPB_True.csv"))
    c.preprocessing(compression=True, expansion=False)
    assert 'sos^tnfa=p38' in c.nodes()



def test_split_node():
    c = CNOGraph()
    c.add_reactions(['a=b', 'a=c', 'c=d', 'b=d', 'd^f=e'])
    c._signals = ['a']
    c._stimuli = ['b']
    c._inhibitors = ['c']
    c.expand_and_gates()
    c.split_node('c',['c1', 'c2'])
    c.split_node('f',['f1', 'f2'])
    assert sorted(c.reactions) == ['a=b', 'a=c1', 'a=c2', 'b=d',
            'b^c1=d','b^c2=d', 'c1=d', 'c2=d', 'd^f1=e', 'd^f2=e']

    c = CNOGraph()
    c.add_reaction("A+B=C")
    c2 = c.copy()

    c._stimuli = ['A']
    c._signals = ['B']
    c._inhibitors = ['C']
    # not yet inplace
    c.relabel_nodes({'A':'a', 'B':'b', 'C':'c'})
    c.split_node('a', ['A1', 'A2'])
    c.split_node('b', ['B1', 'B2'])
    c.split_node('c', ['C1', 'C2'])
    c.merge_nodes(['C1', 'C2'], 'C')
    c.merge_nodes(['B1', 'B2'], 'B')
    c.merge_nodes(['A1', 'A2'], 'A')
    assert c == c2


def test_remove_self_loops():
    c = CNOGraph()
    c.add_reaction("a=b")
    c.add_reaction("b=b")
    c.remove_self_loops()
    c.reactions == ['a=b']

def test_plot():
    # empty graph
    c = CNOGraph()
    c.plot()

    # inh/sti and inh/sig
    c = CNOGraph()
    c.add_reaction("A+B=C")
    c._signals = ['C']
    c._stimuli = ['A', 'B']
    c._inhibitors = ['A', 'C']
    c.plot(rank_method=None, show=False) # speeed up
    fh = tempfile.NamedTemporaryFile(suffix='.svg')

    c.plot(rank_method='all', filename=fh.name, show=False) # speeed up
    c.plot(rank_method='cno', filename=fh.name, show=False) # speeed up

    # node attribute
    c = CNOGraph()
    c.add_reaction("A=B")
    c.node['A']['mycolor'] = '#ff1100'
    # set color for only one node
    c.plot(node_attribute='mycolor')

def test_species():
    c = CNOGraph()
    c.add_reaction("A^B=C")
    assert sorted(c.nodes()) == ['A', 'A^B=C', 'B', 'C']
    assert c.species == ['A', 'B', 'C']



class test_CNOGraph(CNOGraph):

    @classmethod
    def setup_class(klass):
        klass.cnograph = CNOGraph(
                cnodata("PKN-ToyPB.sif"), 
                cnodata("MD-ToyPB.csv"))

    def test_sif(self):
        c2 = CNOGraph(self.cnograph.to_sif())
        print(c2)
        assert c2 == self.cnograph

    def test_relabel_nodes(self):
        c = self.cnograph.copy()
        c.expand_and_gates()
        c.relabel_nodes({'ras':'RAS'})

        # check that the data on node are kept
        c = CNOGraph()
        c.add_reaction("!A+B=C")
        c.node['A']['data'] = [1,2]
        c.relabel_nodes({'A': 'AA'})
        assert c.node['AA']['data'] == [1,2]
        # check that data on edge are kept, in particular the link attribute
        assert c.edge['AA']['C']['link'] == '-'
        
    def test_plotting(self):
        c = self.cnograph.copy()
        c.midas.df += .5
        c.plot(show=True, cmap='heat', colorbar=True)
        c.centrality_closeness()
        c.plot(node_attribute="centrality_closeness", show=False) 
        c.preprocessing()
        c.plot()
    
    def test_edge_attribute(self):
        c = self.cnograph.copy()
        for e in c.edges():
            c.edge[e[0]][e[1]]['penwidth'] = random.random()*10
        for e in c.edges():
            c.edge[e[0]][e[1]]['edge'] = random.random()*10
        c.plot(edge_attribute='edge')

    def test_others(self):
        c = self.cnograph.copy()
        c.verbose = 'INFO'
        c.adjacency_matrix()
        c.summary()
        print(c)
        c.get_same_rank()
        c.lookfor("akt")
        c.lookfor("EGFR")

        c = CNOGraph()
        c.add_node("a=")
        c.add_node("=b")

        # cannot assign anything else than a MIDAS or filemame
        try:
            c.midas = 1
            assert False
        except:
            False
    
        c = CNOGraph()
        c.add_reaction("A=")
        c.add_reaction("A=B")
        c.add_reaction("A=B")
        c.add_reaction("A+C=B")
        c.add_reaction("A^C=B")
        c.add_reaction("a^b^c+a=d")
        c.add_reactions(["A=B", "A=C"])

        # add existing reaction wit another link
        c.add_reaction("!A=B")

        c.get_stats()
        c.reset_edge_attributes()

    def test_centrality(self):
        c = self.cnograph.copy()
        c.centrality_betweeness()
        c.centrality_degree()
        c.centrality_eigenvector()
        c.draw()

    def test_preprocessing(self):
        c = self.cnograph.copy()
        # functional test of the preprocessing steps
        c.preprocessing()

        c = self.cnograph.copy()
        c.cutnonc()
        c.compress()
        c.expand_and_gates(2)
        assert len(c) == 17
        c.expand_and_gates(3) # there is no nodes with 3 inputs so length is the same
        assert len(c) == 17 

        # remove and nodes could be call several times without effect.
        c.remove_and_gates()
        assert len(c) == 12

    def test_swap_edges(self):
        c = self.cnograph.copy()
        c.swap_edges(nswap=10)

    def test_export(self):
        c = self.cnograph.copy()
        fh = tempfile.NamedTemporaryFile(delete=False)
        c.to_sif(fh.name)
        c.to_sbmlqual(fh.name)
        c.to_gexf(fh.name)
        c.to_json(fh.name)
        c.read_json(fh.name)
        fh.delete = True
        fh.close()




def test_link():
    try:
        Link('~')
        assert False
    except:
        assert True
