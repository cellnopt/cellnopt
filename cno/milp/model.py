__author__ = 'ltobalina'

import pulp
from numpy import isnan
from cno.io import Reaction
from itertools import izip


class MILPTrain(object):
    """MILP model for training boolean signalling networks

    For details about the approach see:
    Mitsos A, Melas IN, Siminelakis P, Chairakaki AD, Saez-Rodriguez J, Alexopoulos LG. (2009) Identifying Drug Effects
    via Pathway Alterations using an Integer Linear Programming Optimization Formulation on Phosphoproteomic Data.
    PLoS Comput Biol 5(12): e1000591. doi:10.1371/journal.pcbi.1000591
    """

    def __init__(self, pkn, midas):
        self.pkn = pkn
        self.midas = midas

        self.rxn_raw, self.rxn, self.node, self.experiment = self.initialize_indexing_variables()
        self.R, self.I, self.P = self.initialize_grouping_variables()

        self.model = pulp.LpProblem(name="BoolNetTrain")
        self.y_var, self.z_var, self.x_var = self.define_decision_variables()

    def train(self):
        self.add_problem_constraints()

        # Optimize for error
        error_obj = self.error_objective_expression()
        self.model += error_obj
        sol = self.model.solve()
        if sol is pulp.LpStatusInfeasible:
            raise Exception("Infeasible model.")

        # Optimize for size
        best_fit_value = pulp.value(self.model.objective)
        self.model += error_obj == best_fit_value, 'fit_constraint'
        size_obj = self.size_objective_expression()
        self.model += size_obj
        sol = self.model.solve()
        if sol is pulp.LpStatusInfeasible:
            raise Exception("Infeasible model.")

    def initialize_indexing_variables(self):
        # helper indexing variables
        rxn_raw = self.pkn.reactions
        # replace special characters in reaction names to avoid name interpretation problems in the MILP formulation
        rxn = [rxn.replace("^", "_and_").replace("!", "not_").replace("=", "_eq_") for rxn in rxn_raw]
        node = self.pkn.species
        experiment = [k for k in range(self.midas.nExps)]
        return rxn_raw, rxn, node, experiment

    def initialize_grouping_variables(self):
        # helper group variables
        g_R = dict()  # \mathbf{R}_i: signaling molecules (reactants) for reaction i
        g_I = dict()  # \mathbf{I}_i: inhibitors for reaction i
        g_P = dict()  # \mathbf{P}_i: products for reaction i
        for i, i_raw in izip(self.rxn, self.rxn_raw):
            current_rxn = Reaction(i_raw)

            node_origin = current_rxn.get_signed_lhs_species()
            if node_origin['+']:
                g_R[i] = node_origin['+']
            if node_origin['-']:
                g_I[i] = node_origin['-']

            node_product = current_rxn.rhs
            g_P[i] = node_product

        return g_R, g_I, g_P

    def define_decision_variables(self):
        # Variable declaration
        # y_i: 1 if reaction i is present, 0 otherwise
        # z_i^k: 1 if reaction i takes place in experiment k, 0 otherwise
        # x_j^k: 1 if species j is active in experiment k, 0 otherwise
        y_var = pulp.LpVariable.dicts("y", self.rxn, lowBound=0, upBound=1, cat=pulp.LpBinary)
        z_var = pulp.LpVariable.dicts("z", (self.rxn, self.experiment), lowBound=0, upBound=1, cat=pulp.LpBinary)
        x_var = pulp.LpVariable.dicts("x", (self.node, self.experiment), lowBound=0, upBound=1, cat=pulp.LpBinary)
        return y_var, z_var, x_var

    def add_problem_constraints(self):
        # Nodes under stimuli are set according to the experiments
        # x_j^k = 1, \qquad k=1,\dots,n_e \quad j \in \mathbf{M}^{k,1}
        introduced = self.midas.stimuli
        for k in self.experiment:
            for j in self.midas.names_stimuli:
                value = introduced.get_value(index=introduced.index[k], col=j)
                self.x_var[j][k].bounds(low=value, up=value)

        # Nodes under inhibitory compounds are set to zero when the inhibitor is present
        # x_j^k = 0, \qquad k=1,\dots,n_e \quad j \in \mathbf{M}^{k,0}
        excluded = self.midas.inhibitors
        for k in self.experiment:
            for j in self.midas.names_inhibitors:
                value = excluded.get_value(index=excluded.index[k], col=j)
                if value == 1:
                    self.x_var[j][k].bounds(low=0, up=0)

        # A reaction can only take place if it is possible
        # z_i^k \leq y_i, \qquad i=1,\dots,n_r \quad k=1,\dots,n_e
        for i in self.rxn:
            for k in self.experiment:
                self.model += self.z_var[i][k] <= self.y_var[i]

        # A reaction can only take place if all reagents and no inhibitors are present.
        # z_i^k \leq x_j^k, \qquad i=1,\dots,n_r \quad k=1,\dots,n_e \quad j \in \mathbf{R}_i
        for i in self.rxn:
            for k in self.experiment:
                if i in self.R:
                    for j in self.R[i]:
                        self.model += self.z_var[i][k] <= self.x_var[j][k]
        # z_i^k \leq 1 - x_j^k, \qquad i=1,\dots,n_r \quad k=1,\dots,n_e \quad j \in \mathbf{I}_i
        for i in self.rxn:
            for k in self.experiment:
                if i in self.I:
                    for j in self.I[i]:
                        self.model += self.z_var[i][k] <= 1 - self.x_var[j][k]

        # If a reaction is possible, all reagents are present and no inhibitors are present.
        # z_i^k \geq y_i + \sum_{j \in \mathbf{R}_i} (x_j^k - 1) - \sum_{j \in \mathbf{I_i}} (x_j^k),
        #   \qquad i=1,\dots,n_r \quad k=1,\dots,n_e
        # note: when adding terms to a constraint expression,
        # variables are added to the left of the inequality
        # and constant terms are moved to the right
        # that is why we do -= for \mathbf{R}_i and += for \mathbf{I}_i
        for i in self.rxn:
            for k in self.experiment:
                constraint = self.z_var[i][k] >= self.y_var[i]
                if i in self.R:
                    constraint -= pulp.lpSum((self.x_var[j][k] - 1) for j in self.R[i])
                if i in self.I:
                    constraint += pulp.lpSum(self.x_var[j][k] for j in self.I[i])
                self.model += constraint

        # A species will be formed if some reaction in which it is a product occurs.
        # x_j^k \geq  z_i^k, \qquad i=1,\dots,n_r \quad k=1,\dots,n_e \quad j \in \mathbf{P}_i
        for i in self.rxn:
            for k in self.experiment:
                if i in self.P:
                    j = self.P[i]
                    is_inhibited = False
                    if j in self.midas.names_inhibitors:
                        row_index = self.midas.inhibitors.index[k]
                        is_inhibited = self.midas.inhibitors.get_value(index=row_index, col=j) == 1
                    # add constraint only if the product node has not been manually set
                    if not is_inhibited:
                        self.model += self.x_var[j][k] >= self.z_var[i][k]

        # A species will not be present if all reactions in which it appears as a product do not occur.
        # Note: manipulated species are not considered as products in reactions.
        # x_j^k \leq \sum_{i=1,\dots,n_r; j \in \mathbf{P}_i} z_i^k, \qquad i=1,\dots,n_r \quad k=1,\dots,n_e
        for j in self.node:
            for k in self.experiment:
                # add constraint only if the product node has not been manually set
                is_inhibited = False
                if j in self.midas.names_inhibitors:
                    row_index = self.midas.inhibitors.index[k]
                    is_inhibited = self.midas.inhibitors.get_value(index=row_index, col=j) == 1
                if not is_inhibited:
                    rhs = sum(self.z_var[i][k] for i in self.rxn if i in self.P and j in self.P[i])
                    if rhs is not 0:
                        self.model += self.x_var[j][k] <= rhs

    def error_objective_expression(self):
        # \sum_{k=1,\dots,n_e} \sum_{j \in \mathbf{M}^{k,2}} \alpha_j^k (x_j^{k,m} + (1 - 2 x_j^{k,m}) x_j^k))

        time_start = self.midas.times[0]
        time_end = self.midas.times[1]
        measured = self.midas.df
        measured_start = measured.query("time=="+str(time_start))
        measured_end = measured.query("time=="+str(time_end))

        error_obj = None  # error objective expression initialization
        for k in self.experiment:
            for j in self.midas.names_species:
                value_start = measured_start.get_value(index=measured_start.index[k], col=j)
                value_end = measured_end.get_value(index=measured_end.index[k], col=j)
                value = value_end - value_start
                if not isnan(value):
                    if value < 0:
                        value += 1
                    error_obj += value + (1 - 2*value)*self.x_var[j][k]
        return error_obj

    def size_objective_expression(self):
        # \sum_{i=1,\dots,n_r} \beta_i y_i
        size_obj = pulp.lpSum(self.y_var)
        return size_obj
