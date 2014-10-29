from cno.io import sif2asp
from cno import cnodata
import os

cleanup = "PKN-ToyMMB.net"

def test_sif2aspnet_py():
    filename = "PKN-ToyMMB.sif"
    model = cnodata(filename)
    c = sif2asp.SIF2ASP(model)
    c.to_net(filename.replace(".sif", ".net"))
    os.remove(cleanup) 
    print(c)

