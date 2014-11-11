from cno import cnodata, CNOGraph, SIF
from cno.io.sbmlqual import SBMLQual
from easydev import TempFile

from cno.io.cnograph import CNOGraph


def test_sbmlqual_class():
    s = SBMLQual()
    try:
        s.to_sbmlqual(1)
        assert False
    except:
        assert True
    c = CNOGraph()
    c.add_reaction("A=B")
    s.to_sbmlqual(c)


def test_sif_vs_sbml_datasets():
    for identifier in ['ToyPCB', 'ToyPB', 'ToyPB_True', 'ToyMMB', 
            'LiverDREAM', 'ExtLiverPCB']:
        # should be enough for now
        yield sbmlqual_from_datasets, identifier


def sbmlqual_from_datasets(identifier):

    # a simple model
    s1 = SIF()
    s2 = SIF(cnodata("PKN-" + identifier + ".sif"))
    fh = TempFile()
    s2.to_sbmlqual(fh.name)
    s1.read_sbmlqual(fh.name)
    fh.delete()
    assert s1 == s2
    s3 = SIF(cnodata("PKN-" + identifier + ".xml"))
    assert s1 == s3 and s2 == s3



def test_read_write_from_cnograph():
    c  = CNOGraph(cnodata("PKN-ToyPB.sif"))
    fh = TempFile(suffix='.xml')
    c.to_sbmlqual(fh.name)
    c2 = CNOGraph(fh.name)
    assert c == c2
    fh.delete()

    c  = CNOGraph(cnodata("PKN-ToyPB.sif"))
    c.expand_and_gates()
    fh = TempFile(suffix='.xml')
    c.to_sbmlqual(fh.name)
    c2 = CNOGraph(fh.name)
    fh.delete()
    assert c == c2

def test_simple_sbmlqual():
    # a simple  example with simple OR, simple link, mix of OR and AND and single ANd 
    c = CNOGraph()
    c.add_reaction("!A=C")
    c.add_reaction("C=D")
    c.add_reaction("B=C")
    c.expand_and_gates()
    c.add_reaction("a1=b")
    c.add_reaction("a2=b")
    c.add_reaction("D^b=E")


    fh = TempFile(suffix='.xml')
    c.to_sbmlqual(fh.name)
    c2 = CNOGraph(fh.name)
    fh.delete()
    assert c == c2



def test_toypb_bioservices():
    from bioservices import biomodels
    b = biomodels.BioModels()

    sbml = b.getModelSBMLById('MODEL1305240000')

    fh =TempFile(suffix='.xml')
    with open(fh.name, 'w') as fh:
        fh.write(sbml.encode('utf-8'))

    c = CNOGraph(fh.name)

    c2 = CNOGraph(cnodata("PKN-ToyPB.sif"))
    assert c == c2



