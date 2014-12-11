import numpy as np
import pylab
import pandas as pd
from biokit.rtools import RSession
from cno.core import CNORBase
from easydev import TempFile
__all__ = ['GTTBool']


class GTTBool(CNORBase):
    """


    ::

        from cno import *
        c = cnorbool.CNORbool(cnodata("PKN-ToyMMB.sif"), cnodata("MD-ToyMMB.csv"), verboseR=False)
        c.optimise(reltol=0.5)
        c.optimise(reltol=0.5)
        g   = gtt.GTTBool(c._model, c.data, c.models, c.results.scores)
        d = g.get_gtt()
    """
    def __init__(self, model, data, models, scores, verboseR=False):
        """
        Note that once, grouped, the scores should be identical albeit the
        model size

        [scores[i] for i in grouped.groups.values()[10]]


        :param model: a instance of :class:`CNOGraph`
        :param data: an instance of :class:`XMIDAS`
        :param models: an instance of compatible :class:`Models`
        :param scores: the scores of each model.





        """
        CNORBase.__init__(self, verboseR)
        self.models = models
        self.scores = scores
        self.model = model
        self.data = data # a MIDAS file


    def _init(self):
        fhmodel = TempFile()
        fhdata = TempFile()
        self.model.to_sif(fhmodel.name)
        self.data.to_midas(fhdata.name)
        self.session.run("library(CellNOptR)")
        self.session.run('model=readSIF("%s")' % fhmodel.name)
        self.session.run('cnolist=CNOlist("%s")' % fhdata.name)

    def _get_sim(self, bs):
        self.session.bs1 = bs
        script = """
        png()
        output = cutAndPlot(cnolist, model, list(bs1), plotPDF=F)
        dev.off()
        """
        self.session.run(script)
        res = self.session.output['simResults'][0]
        res = list(res['t0'].flatten() ) + list(res['t1'].flatten())
        return res

    def get_gtt(self):
        print("init R library")
        self._init()
        N = len(self.models)
        from easydev import progress_bar
        b = progress_bar(N)
        d = {}
        for i in range(0, N):
            res = np.array(self._get_sim(self.models.df.ix[i].values))
            b.animate(i, N)
            d[i] = res

        df = pd.DataFrame(d).transpose()
        grouped = df.groupby(list(df.columns))
        pylab.hist([len(this) for this in grouped.groups.values()], 100)
        return {'simulation': d, 'grouped':grouped}


