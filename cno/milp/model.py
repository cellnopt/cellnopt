# -*- python -*-
#
# This file is part of the CNO package
#
#  Copyright (c) 2014 - EMBL-EBI
#
#  File author(s): Luis Tobalina
#    Thomas Cokelaer (cokelaer@ebi.ac.uk)
#
#  Distributed under the GLPv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: www.cellnopt.org
#
##############################################################################
import sys

import pulp
from numpy import isnan
from cno.io import Reaction
from cno.core.base import CNOBase

from cno.core.params import OptionsBase


class MILPTrain(CNOBase):
    """Mixed Integer Linear Program (MILP) model for training boolean signalling networks.

    This class builds an optimization model from a given prior knowledge
    network (pkn) and data in midas format. The aim is to select a subset of
    interactions from the given pkn that are able to explain the observed data
    under a boolean network framework.

    The problem is solved in two steps. First, the algorithm searches for a
    network with the best possible fit. Second,the algorithm searches for
    the smallest possible network that returns the best fit found.

    Example::

        >>> from cno import cnodata, CNOGraph, XMIDAS
        >>> import cno.milp.model

        >>> filename_pkn = cnodata("PKN-ToyMMB.sif")
        >>> filename_midas = cnodata("MD-ToyMMB.csv")
        >>> pkn = CNOGraph(filename_pkn, filename_midas)
        >>> midas = pkn.midas

        >>> pkn.compress()
        >>> pkn.expand_and_gates()

        >>> model = cno.milp.model.MILPTrain(pkn, midas)
        >>> model.train()
        >>> #model.get_rxn_solution()

    The problem is implemented using PuLP, a linear programming toolkit
    for python than can interface with different linear programming solvers,
    both commercial (e.g. CPLEX, Gurobi) and open-source (e.g. CBC).
    The pulp.LpProblem() instance is stored in MILPTrain.model attribute.
    The user can select a different solver by using the setSolver()
    method provided in the pulp.LpProblem() instance. For more details,
    see the PuLP documentation.

    For details about the approach see:

    * Mitsos A, Melas IN, Siminelakis P, Chairakaki AD, Saez-Rodriguez J,
      Alexopoulos LG. (2009) Identifying Drug Effects via Pathway Alterations
      using an Integer Linear Programming Optimization Formulation on
      Phosphoproteomic Data. PLoS Comput Biol 5(12): e1000591.
      doi:10.1371/journal.pcbi.1000591

    .. note:: the constraint explained in equation number 2 in the paper is
        not included in this implementation. This constraint would allow the
        modeller to limit the combinations of connectivities considered. It
        is thus not essential for solving the main problem.
    """

    def __init__(self, pkn, midas, verbose=False):
        """Initialization function.

        :param CNOGraph pkn: prior knowledge network. see :class:`~cno.io.cnograph.CNOGraph`
        :param XMIDAS midas: experimental data. See :class:`~cno.io.midas.XMIDAS`
        """
        super(MILPTrain, self).__init__(pkn, midas, verbose=verbose)

        self.rxn_raw, self.rxn, self.node, self.experiment = self.initialize_indexing_variables()
        self.R, self.I, self.P = self.initialize_grouping_variables()

        self.model = pulp.LpProblem(name="BoolNetTrain")
        self.y_var, self.z_var, self.x_var = self.define_decision_variables()

    def _get_pknmodel(self):
        return self._pknmodel

    def _set_pknmodel(self, new_pkn):
        self._pknmodel = new_pkn
        self._reset_class_attributes()
        self._reset_problem()

    pknmodel = property(_get_pknmodel, _set_pknmodel,
                        doc="getter/setter to the model")

    def _get_midas(self):
        return self._data

    def _set_midas(self, new_midas):
        self._data = new_midas
        self._reset_class_attributes()
        self._reset_problem()

    midas = property(_get_midas, _set_midas, doc="getter/setter to the data")

    def change_pkn_and_midas(self, new_pkn, new_midas):
        """Change information about the stored pkn and midas.

        This function changes both attributes, the pkn and the midas.

        :param CNOGraph new_pkn: new prior knowledge network.
        :param XMIDAS new_midas: new experimental data.
        :return:
        """
        self._pknmodel = new_pkn
        self._midas = new_midas
        self._reset_class_attributes()
        self._reset_problem()

    def _reset_class_attributes(self):
        r"""Reset class attributes.

        This function recalculates rxn_raw, rxn, node, experiment, R, I, P,
        :math:`y_var`, :math:`z_var` and :math:`x_var`.
        It is designed to be called when either pkn or midas are
        changed by the user.
        """
        self.rxn_raw, self.rxn, self.node, self.experiment = self.initialize_indexing_variables()
        self.R, self.I, self.P = self.initialize_grouping_variables()

        self.y_var, self.z_var, self.x_var = self.define_decision_variables()

    def _reset_problem(self):
        """Reset the optimization problem to an empty problem.

        This function removes all the variables and constraints, as well as
        the objective function, from the  pulp.LpProblem() instance.
        """
        # remove objective function
        self.model.objective = None
        # remove constraints
        self.model.constraints.clear()
        # remove variables
        # Note: I cannot find any method in the pulp.LpProblem() class to
        # remove the attached variables. This two lines
        # make the trick, but may not be optimal.
        self.model._variable_ids.clear()
        self.model._variables = []

    def optimise(self):
        return self.train()

    def train(self):
        """Initialize and solve the optimization problem.

        The problem is solved in two steps. First, the algorithm searches for
        a network with the best possible fit. Second, the algorithm searches
        for the smallest possible network that returns the best fit found.
        """
        # reset constraint set
        self.model.constraints.clear()

        # add problem constraints
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
        """Initialize indexing variables used in the formulation of the model."""
        # helper indexing variables
        rxn_raw = self.pknmodel.reactions
        # replace special characters in reaction names to avoid name interpretation problems in the MILP formulation
        rxn = [rxn.replace("^", "_and_").replace("!", "not_").replace("=", "_eq_") for rxn in rxn_raw]
        node = self.pknmodel.species
        experiment = [k for k in range(self.midas.nExps)]
        return rxn_raw, rxn, node, experiment

    def initialize_grouping_variables(self):
        """Initialize grouping variables used in the formulation of the model.

        * :math:`\mathbf{R}_i`: signaling molecules (reactants) for reaction i.
        * :math:`\mathbf{I}_i`: inhibitors for reaction i.
        * :math:`\mathbf{P}_i`: products for reaction i.
        """
        # helper group variables
        g_R = dict()  # \mathbf{R}_i: signaling molecules (reactants) for reaction i
        g_I = dict()  # \mathbf{I}_i: inhibitors for reaction i
        g_P = dict()  # \mathbf{P}_i: products for reaction i
        for i, i_raw in zip(self.rxn, self.rxn_raw):
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
        """Decision variables used in the formulation of the model.

        * :math:`y_i`: 1 if reaction i is present, 0 otherwise.
        * :math:`z_i^k`: 1 if reaction i takes place in experiment k, 0 otherwise.
        * :math:`x_j^k`: 1 if species j is active in experiment k, 0 otherwise.
        """
        # Variable declaration
        y_var = pulp.LpVariable.dicts("y", self.rxn, lowBound=0, upBound=1, cat=pulp.LpBinary)
        z_var = pulp.LpVariable.dicts("z", (self.rxn, self.experiment), lowBound=0, upBound=1, cat=pulp.LpBinary)
        x_var = pulp.LpVariable.dicts("x", (self.node, self.experiment), lowBound=0, upBound=1, cat=pulp.LpBinary)
        return y_var, z_var, x_var

    def add_problem_constraints(self):
        """Add modelling constraints to the optimization problem.

        This function adds constraints described in equations (3) to (10) of
        the paper Mitsos A, Melas IN, Siminelakis P, Chairakaki AD,
        Saez-Rodriguez J, Alexopoulos LG. (2009) Identifying Drug Effects via
        Pathway Alterations using an Integer Linear Programming Optimization
        Formulation on Phosphoproteomic
        Data. PLoS Comput Biol 5(12): e1000591. doi:10.1371/journal.pcbi.1000591
        """

        # equation number 10 in Mitsos et al. 2009
        # Nodes under stimuli are set according to the experiments
        # x_j^k = 1, \qquad k=1,\dots,n_e \quad j \in \mathbf{M}^{k,1}
        introduced = self.midas.stimuli
        for k in self.experiment:
            for j in self.midas.names_stimuli:
                value = introduced.get_value(index=introduced.index[k], col=j)
                self.x_var[j][k].bounds(low=value, up=value)

        # equation number 9 in Mitsos et al. 2009
        # Nodes under inhibitory compounds are set to zero when the inhibitor is present
        # x_j^k = 0, \qquad k=1,\dots,n_e \quad j \in \mathbf{M}^{k,0}
        excluded = self.midas.inhibitors
        for k in self.experiment:
            for j in self.midas.names_inhibitors:
                value = excluded.get_value(index=excluded.index[k], col=j)
                if value == 1:
                    self.x_var[j][k].bounds(low=0, up=0)

        # equation number 3 in Mitsos et al. 2009
        # A reaction can only take place if it is possible
        # z_i^k \leq y_i, \qquad i=1,\dots,n_r \quad k=1,\dots,n_e
        for i in self.rxn:
            for k in self.experiment:
                self.model += self.z_var[i][k] <= self.y_var[i], 'C3_{}_{}'.format(i, k)

        # equation number 4 in Mitsos et al. 2009
        # A reaction can only take place if all reagents and no inhibitors are present.
        # z_i^k \leq x_j^k, \qquad i=1,\dots,n_r \quad k=1,\dots,n_e \quad j \in \mathbf{R}_i
        for i in self.rxn:
            for k in self.experiment:
                if i in self.R:
                    for j in self.R[i]:
                        self.model += self.z_var[i][k] <= self.x_var[j][k], 'C4_{}_{}_{}'.format(i, k, j)
        # equation number 5 in Mitsos et al. 2009
        # z_i^k \leq 1 - x_j^k, \qquad i=1,\dots,n_r \quad k=1,\dots,n_e \quad j \in \mathbf{I}_i
        for i in self.rxn:
            for k in self.experiment:
                if i in self.I:
                    for j in self.I[i]:
                        self.model += self.z_var[i][k] <= 1 - self.x_var[j][k], 'C5_{}_{}_{}'.format(i, k, j)

        # equation number 6 in Mitsos et al. 2009
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
                self.model += constraint, 'C6_{}_{}'.format(i, k)

        # equation number 7 in Mitsos et al. 2009
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
                        self.model += self.x_var[j][k] >= self.z_var[i][k], 'C7_{}_{}_{}'.format(i, k, j)

        # equation number 8 in Mitsos et al. 2009
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
                        self.model += self.x_var[j][k] <= rhs, 'C8_{}_{}'.format(j, k)

    def error_objective_expression(self):
        r"""Define the error objective function.

        .. math::

            \sum_{k=1,\dots,n_e} \sum_{j \in \mathbf{M}^{k,2}} \alpha_j^k (x_j^{k,m} + (1 - 2 x_j^{k,m}) x_j^k))

        :return: pulp.LpAffineExpression

        """
        # equation number 1 (first part) in Mitsos et al. 2009
        time_start = self.midas.times[0]
        time_end = self.midas.times[1]
        measured = self.midas.df
        measured_start = measured.query("time==" + str(time_start))
        measured_end = measured.query("time==" + str(time_end))

        error_obj = None  # error objective expression initialization
        for k in self.experiment:
            for j in self.midas.names_species:
                value_start = measured_start.get_value(index=measured_start.index[k], col=j)
                value_end = measured_end.get_value(index=measured_end.index[k], col=j)
                value = value_end - value_start
                if not isnan(value):
                    if value < 0:
                        value += 1
                    error_obj += value + (1 - 2 * value) * self.x_var[j][k]
        return error_obj

    def size_objective_expression(self):
        r"""Define the network size objective function.

        .. math::
            \sum_{i=1,\dots,n_r} \beta_i y_i
        """
        # equation number 1 (second part) in Mitsos et al. 2009
        size_obj = pulp.lpSum(self.y_var)
        return size_obj

    def get_rxn_solution(self):
        """Get which reactions are part of the optimized network.

        :return: dictionary with reaction names as keys and 1 or 0 as value i
            indicating if the reaction is present in the optimized network or not.
        """
        rxn_sol = dict.fromkeys(self.rxn_raw, 0)
        for i, i_raw in zip(self.rxn, self.rxn_raw):
            rxn_sol[i_raw] = self.y_var[i].value()
        return rxn_sol



def standalone(args=None):
    """This function is used by the standalone application called cellnopt_boolean

    ::

        cno_milp --help

    """
    if args is None:
        args = sys.argv[:]

    user_options = OptionsMILP()

    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    else:
        options = user_options.parse_args(args[1:])

    if options.onweb is True or options.report is True:
        o = MILPTrain(options.pknmodel, options.data, verbose=options.verbose
             )

    if options.onweb is True:
        o.optimise()
        o.onweb()
    elif options.report is True:
        o.optimise()
        o.report()
    else:
        from easydev.console import red
        print(red("No report requested; nothing will be saved or shown"))
        print("use --on-web or --report options")


class OptionsMILP(OptionsBase):
    def __init__(self):
        prog = "cno_milp"
        version = prog + " v1.0 (Thomas Cokelaer @2014)"
        super(OptionsMILP, self).__init__(version=version, prog=prog)

if __name__ == "__main__":
    standalone()
