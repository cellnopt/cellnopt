"""Example of cross validation that could be used on 
steady case.

"""

from cno.boolean.steady import Steady







class CrossValidation(object):
    def __init__(self, pknmodel, midas):

        self.steady = Steady(pknmodel, midas)
        print('Perform one optimisation to have a good guess')
        self.steady.preprocessing()
        self.steady.optimise(verbose=True)
        self.best_bitstring = self.steady.best_bitstring
        print(self.steady.results.scores[-1])
        self.experiments = self.steady.midas.experiments
        self.midas = midas
        self.pknmodel = pknmodel


    def one_fold(self):
        # for each experiment, leave it out, optimise and get errors
        # Then, keep the optimised reactions

        self.scores_one_out = []
        self.bs_one_out = []
        self.scores_test = []
        for experiment in self.experiments.index:
            steady = Steady(self.pknmodel, self.midas)
            steady.midas.remove_experiments([experiment])
            steady._init_model_data()
            steady.preprocessing()
            steady.optimise(guess=self.best_bitstring)
            print experiment, steady.results.scores[-1]
            self.scores_one_out.append(steady.results.scores[-1])
            self.bs_one_out.append(steady.best_bitstring[:])

            # now compute errors of the test experiement
            # keeping only the test experiment, simulate and compute score
            steady_test = Steady(self.pknmodel, self.midas)
            steady_test.preprocessing()
            steady_test.midas.remove_experiments([exp for exp in self.experiments.index if exp!=experiment])
            steady_test._init_model_data()
            steady_test.simulate(reactions=steady_test.parameters2reactions(steady.best_bitstring[:]))
            self.scores_test.append(steady_test.score())
            print steady_test.score()
            print("----")
            self.steady_test = steady_test







