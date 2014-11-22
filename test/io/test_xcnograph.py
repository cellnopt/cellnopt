from cno.io import CNOGraph, XMIDAS, SIF
from cno import cnodata, getdata, CNOError
import tempfile
import random
from cno.io import XCNOGraph


sif = getdata("PKN-test_cnograph.sif")
midas = getdata("MD-test.csv")




def test_xcnograph():


    c = XCNOGraph(sif, midas)
    c.plot_adjacency_matrix()
    c.plot_feedback_loops_species()
    c.plot_feedback_loops_histogram()
    c.degree_histogram(show=True)
    c.hcluster()
    c.plot_degree_rank(loc='upper right', layout='circular')
    c.plot_degree_rank(loc='other', layout='spring')
    c.plot_in_out_degrees()
    c.dependency_matrix()

    c = XCNOGraph()
    c.plot_feedback_loops_species()


def test_poisson():
    c = XCNOGraph()
    c.random_poisson_graph(n=100, mu=2.5, remove_unconnected=False)
    # test also recursive compression
    c.expand_and_gates()
    c.compress()

