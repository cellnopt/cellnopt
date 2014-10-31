from cno.io import CNA
from cno.testing import getdata


import numpy
import tempfile
import os

from easydev import gsf


def test_cna_reaction():
   a = CNA(getdata("test_cna_reactions"), verbose=True)
   assert len(a.reactions) == 123
   assert len(a.species) == 94

   f = tempfile.NamedTemporaryFile()
   a.to_sif(f.name)


   f = tempfile.NamedTemporaryFile()
   a.to_sif(f.name)


def test1():
    r = CNA()
    r.add_reaction("a+b=g")
    r.add_reaction("a+b=c")
    assert r.species ==  ['a', 'b', 'c', 'g']


    r.to_sif("test.sif")  # this creates the nodes attributes that change behabour of specID
    r.remove_reaction("a+b=c")
    assert r.species ==  ['a', 'b', 'g']
    os.remove("test.sif")
