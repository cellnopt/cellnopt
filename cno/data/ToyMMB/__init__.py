from cellnopt.core import XMIDAS, CNOGraph

pknmodel = CNOGraph("PKN-ToyMMB.sif")
data = XMIDAS("MD-ToyMMB.csv")
description = open("README.rst").read()
