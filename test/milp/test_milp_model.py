import cno.milp.model
import pulp

from cno import cnodata
from cno.io import CNOGraph, XMIDAS


class TestMILPTrain(object):
    def setup(self):
            filename_pkn = cnodata("PKN-ToyMMB.sif")
            filename_midas = cnodata("MD-ToyMMB.csv")

            self.pkn = CNOGraph(filename_pkn, filename_midas)
            self.midas = self.pkn.midas

    def teardown(self):
        pass

    def test_init(self):
        model = cno.milp.model.MILPTrain(self.pkn, self.midas)
        nreactions = len(self.pkn.reactions)
        nspecies = len(self.pkn.species)
        nexperiments = self.midas.nExps
        assert len(model.y_var) == nreactions
        assert len(model.z_var) == nreactions
        assert len(model.z_var[model.rxn[0]]) == nexperiments
        assert len(model.x_var) == nspecies
        assert len(model.x_var[model.node[0]]) == nexperiments

    def test_expanded_ToyMMB(self):
        self.pkn.compress()
        self.pkn.expand_and_gates()

        model = cno.milp.model.MILPTrain(self.pkn, self.midas)
        model.train()

        eps = 1e-12
        model_error = pulp.value(model.error_objective_expression())
        model_size = pulp.value(model.size_objective_expression())

        assert abs(model_error-10.02) < eps
        assert abs(model_size-9) < eps
