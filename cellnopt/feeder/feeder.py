# -*- python -*-
#
#  This file is part of cellnopt software
#
#  Copyright (c) 2014 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: http://github.com/cellnopt
#
##############################################################################
from biokit import rtools


class Feeder(object):
    def __init__(self):
        self.Rscript = """
            library(CNORfeeder)
            pknmodel = readSIF("%(pkn)s")
            cnolist = CNOlist("%(data)s")

            model = preprocessing(cnolist, pknmodel,
                                  compression=%(compression)s,
                                  expansion=%(expansion)s,
                                  maxInputsPerGate=%(maxInputsPerGate)s)
            # find new links
            BTable = makeBTables(CNOlist=cnolist, k=%(k)s,
                                measErr=c(%(err1)s,%(err2)s))
            alllinks = linksRanking(CNOlist=cnolist,
                                    measErr=c(%(err1)s,%(err2)s))
            modelIntegr = mapBTables2model(BTable=BTable, model=model,
                                           allInter=F)
            newlinks = modelIntegr$reacID[modelIntegr$indexIntegr]
            """
        self.rsession = rtools.RSession()

    def run(self, model, data, k=2, compression=True,
            expansion=True, maxInputsPerGate=3, err1=0.1, err2=0,
            verbose=True):

        cmd = self.Rscript % {'pkn': model,
                              'data': data,
                              'k': k,
                              'compression': rtools.bool2R(compression),
                              'expansion': rtools.bool2R(expansion),
                              'err1': err1,
                              'err2': err2,
                              'maxInputsPerGate': maxInputsPerGate}

        if verbose:
            print(cmd)
        self.rsession.run(cmd)
        self.alllinks = self.rsession.get('alllinks')
        self.newlinks = self.rsession.get('newlinks')

    def summary(self):
        print("Found {0} new links".format(len(self.newlinks)))
        print(self.newlinks)
