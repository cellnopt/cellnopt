__author__ = 'ltobalina'

import pulp
from numpy import isnan


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
        rxn_raw = [e[0]+"->"+e[1] for e in self.pkn.edges_iter() if "=" not in e[1]]
        rxn = [rxn.replace("^", "_and_").replace("!", "not_").replace("=", "_e_") for rxn in rxn_raw]
        node = [n for n in self.pkn.nodes_iter() if "=" not in n]
        experiment = [k for k in range(self.midas.nExps)]
        return rxn_raw, rxn, node, experiment

    def initialize_grouping_variables(self):
        # helper group variables
        g_R = dict()  # R_i: signaling molecules (reactants) for reaction i
        g_I = dict()  # I_i: inhibitors for reaction i
        g_P = dict()  # P_i: products for reaction i
        for i in self.rxn_raw:
            node_origin, node_product = i.split("->")
            g_P[i] = node_product
            if "=" not in node_origin:
                link_type = self.pkn.get_edge_data(node_origin, node_product)['link']
                if link_type == '+':
                    g_R[i] = [node_origin]
                elif link_type == '-':
                    g_I[i] = [node_origin]
            else:
                for pred in self.pkn.predecessors(node_origin):
                    link_type = self.pkn.get_edge_data(pred, node_origin)['link']
                    if link_type == '+':
                        if i not in g_R:
                            g_R[i] = []
                        g_R[i].append(pred)
                    elif link_type == '-':
                        if i not in g_I:
                            g_I[i] = []
                        g_I[i].append(pred)

        return g_R, g_I, g_P

    def define_decision_variables(self):
        # Variable declaration
        # y_i: 1 if reaction i present, 0 otherwise
        # z_i^k: 1 if reaction i takes place in experiment k, 0 otherwise
        # x_j^k: 1 if species j is active in experiment k, 0 otherwise
        y_var = pulp.LpVariable.dicts("rxn", self.rxn, lowBound=0, upBound=1, cat=pulp.LpBinary)
        z_var = pulp.LpVariable.dicts("rxn", (self.rxn, self.experiment), lowBound=0, upBound=1, cat=pulp.LpBinary)
        x_var = pulp.LpVariable.dicts("n", (self.node, self.experiment), lowBound=0, upBound=1, cat=pulp.LpBinary)
        return y_var, z_var, x_var

    def add_problem_constraints(self):
        # Nodes under stimuli are set according to the experiments
        # x_j^k = 1    k=1,...,n_e    j in M^(k,1)
        introduced = self.midas.stimuli
        for k in self.experiment:
            for j in self.midas.names_stimuli:
                value = introduced.get_value(index=introduced.index[k], col=j)
                self.x_var[j][k].bounds(low=value, up=value)

        # Nodes under inhibitory compounds are set to zero when the inhibitor is present
        # x_j^k = 0    k=1,...,n_e    j in M^(k,0)
        excluded = self.midas.inhibitors
        for k in self.experiment:
            for j in self.midas.names_inhibitors:
                value = excluded.get_value(index=excluded.index[k], col=j)
                if value == 1:
                    self.x_var[j][k].bounds(low=0, up=0)

        # A reaction can only take place if it is possible
        # z_i^k <= y_i    i=1,...,n_r    k=1,...,n_e
        for i in self.rxn:
            for k in self.experiment:
                self.model += self.z_var[i][k] <= self.y_var[i]

        # A reaction can only take place if all reagents and no inhibitors are present.
        # z_i^k <= x_j^k    i=1,...,n_r    k=1,...,n_e    j in R_i
        for i in self.rxn:
            for k in self.experiment:
                if i in self.R:
                    #j = self.R[i]
                    #self.model += self.z_var[i][k] <= self.x_var[j][k]
                    for j in self.R[i]:
                        self.model += self.z_var[i][k] <= self.x_var[j][k]
        # z_i^k <= 1 - x_j^k    i=1,...,n_r    k=1,...,n_e    j in I_i
        for i in self.rxn:
            for k in self.experiment:
                if i in self.I:
                    #j = self.I[i]
                    #self.model += self.z_var[i][k] <= 1 - self.x_var[j][k]
                    for j in self.I[i]:
                        self.model += self.z_var[i][k] <= 1 - self.x_var[j][k]

        # If a reaction is possible, all reagents are present and no inhibitors are present.
        # z_i^k >= y_i + sum(x_j^k - 1)_(j in R_i) - sum(x_j^k)_(j in I_i)    i=1,...,n_r    k=1,...,ne
        # note: when adding terms to a constraint expression,
        # variables are added to the left of the inequality
        # and constant terms are moved to the right
        # that is why we do -= for R_i and += for I_i
        for i in self.rxn:
            for k in self.experiment:
                constraint = self.z_var[i][k] >= self.y_var[i]
                if i in self.R:
                    #j = self.R[i]
                    #constraint -= self.x_var[j][k]-1
                    constraint -= pulp.lpSum(self.x_var[j][k] for j in self.R[i]) - 1
                if i in self.I:
                    #j = self.I[i]
                    #constraint += self.x_var[j][k]
                    constraint += pulp.lpSum(self.x_var[j][k] for j in self.I[i])
                self.model += constraint

        # A species will be formed if some reaction in which it is a product occurs.
        # x_j^k >= z_i^k    i=1,...,n_r    k=1,...,n_e    j in P_i
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
        # x_j^k <= sum(z_i^k)_(i=1,...,n_r; j in P_i)    i=1,...,n_r    k=1,...,n_e
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
        # sum(sum(alpha_j^k * (x_j^(k,m) + (1 - 2*x_j^(k,m))*x_j^k))_(j in M^(k,2)))_(k=1,...,n_e)

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
        # sum(beta_i*y_i)_(i=1,...,n_r)
        size_obj = pulp.lpSum(self.y_var)
        return size_obj
