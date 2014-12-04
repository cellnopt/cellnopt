# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 09:09:46 2014

@author: cokelaer
"""


class Fitness(object):
    """Not to be used


    """
    def __init__(self, midas, bitstring, reacID, sizeFac=1e-4, NAFac=1):
        """midas must contain a sim attribute


        """
        self.midas = midas
        self.sim = midas.sim
        self.bitstring = bitstring
        self.reacID = reacID

        self.sizeFac = sizeFac

        self.NAFac = NAFac

    def get_score(self, time=30, verbose=True):
        """

        !!! issue with time0

        If we take it into account (deviationPen)
        """
        # number of NA in the data at time 30
        assert time in self.midas.times
        na_indata = self.midas.df.query("time==@time").isnull().sum().sum()
        print(na_indata)

        # not clear in cellnopt if t0 is taken into account. It looks like it is not
        # number of data points
        nDataPts = (self.midas.df.shape[0] * self.midas.df.shape[1])
        nData = float(nDataPts - na_indata) # number of valid data points for normalisation


        reactions = [reac for bit,reac in zip(self.bitstring, self.reacID)]
        inputs = [reac.split("=")[0] for reac in reactions]
        nInTot = sum([len(this.split('+')) for this in inputs])



        reactions = [reac for bit,reac in zip(self.bitstring, self.reacID) if bit==1]
        inputs = [reac.split("=")[0] for reac in reactions]
        nInputs = sum([len(this.split('+')) for this in inputs])

        deviationPen = ((self.midas.sim - self.midas.df)**2).sum().sum() # if t0 taken into account divide by 2.
        deviationPen/=nData


        # artefact of CNOR code where multiply by N at time 0 (ignoring NA)
        # and dividing by N at ti removing NA . Why ?
        sizePen = (nDataPts * self.sizeFac * nInputs)/nInTot
        sizePen /= nData


        # unnormalised NA penalty from simulation
        NAPen = self.NAFac * self.midas.sim.isnull().sum().sum()
        NAPen /= nData
        score = deviationPen + NAPen + sizePen

        print("SCore: {} (deviation={} - NA={} - size={})".format(score,
              deviationPen*nData, NAPen*nData, sizePen*nData))
        print("SCore: {} (deviation={} - NA={} - size={})".format(score, deviationPen, NAPen, sizePen))

        # here residual is computed across all data and normalised accross all data
        N = float(self.midas.df.values.size)
        # FIXME ignore NA ?
        print("Residual errors = {}".format(self.midas.get_residual_errors().sum().sum()/N    ))

        print("nData = {}".format(nData))
        print("nInputs = {}".format(nInputs))
        print("nInTot = {}".format(nInTot))
        print("na_indata = {}".format(na_indata))

        # they could be slight differences with cell not when t0 is taken into account
        #

        return (score)

