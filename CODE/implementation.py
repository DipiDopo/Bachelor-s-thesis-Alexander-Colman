# Author: Alexander Colman, student Physics at the University of Antwerp
# Contact: Alexander dot colman1 AT gmail dot com
# Latest update: 10/05/2025

from datetime import datetime
from gurobipy import Model, GRB, quicksum, Column
import igraph as ig
import matplotlib.pyplot as plt
import random


# This program is an implementation of the branch-and-bound algorithm, using
# the Simplex method and column generation.

# It is designed to solve a linear programming (LP) problem with a focus on
# path generation and branching. There are multiple ways to run this code;
# the initial focal point of this code was to implement the example from the
# paper "A Primer in Column Generation", by J. Desrosiers and M.E. Lubbecke.
# 2005. The functions for running that example can be found in the functions
# below the ColumnGeneration class (BB0, BB1, BB2, BB3 and BB4).

# This code can also be used to solve randomly generated DAGs (directed
# acyclic graphs). Note that the complexity of the LP increases exponentially
# when adding more nodes and arcs.

# The file Alexander_Colman_Bachelorproef_Runner.py can be used to run the
# code concerning the Primer example.

# Finally, some research was done to compare the computation times of
# different branching tactics. The two commented functions at the bottom
# of this code were initially used for this, still using the example from the
# Primer. However, the example is not complex enough; only 6 nodes and 10 arcs
# are incorporated in the LP, which means the computation times are relatively
# small.
# To solve this issue, a new LP is necessary; the code in the file
# Alexander_Colman_Time_Research.py created a model with random costs and
# times between 1 and 20. This model is then used to benchmark the computation
# times for the 3 branching tactics. These tactics are comprised of:
# - Branching on 2 arcs, but inserting them in 2 seperate constraints
# - Inserting 2 arcs into the same constraint
# - Branching on a whole path
# The research compares the computation times for enforcing and omitting (2 of
# the arcs of the) 2 non-integer paths. Additionally, the effect of altering
# the constraint value of enforcing the whole path is explored: initially,
# this is set to 100, but the values 99, 98, 97, 96, 50 and 1 are examined as
# well.


class ColumnGeneration:
    """ A class to calculate the path with the least cost, using the column
   generation algorithm. """

    def __init__(self, arcs: list, times: dict, costs: dict,
                 time_limit: int = 14, dummy_nodes: list = [False, False],
                 initial_cost: int = 1000,
                 objective: str = "MINIMIZE") -> None:
        """ Generates a path with the least cost, while taking all constraints
        into account, using column generation algorithm. The Simplex method is
        used to get an initial solution so the algorithm can start. The
        Branch-and-Bound algorithm is an optional addition to the column
        generation algorithm, to solve the potential "problem" of fractional
        paths.
        ======================================================================
        ATTRIBUTES:
        arcs:
        This is a list of tuples contraining every arc you wish there to be
        between nodes. Please make sure a path is possible from a sink node to
        a source node.
        ----------------------------------------------------------------------
        times:
        This is a dictionary containing the time values (a.k.a. "time-costs")
        per arc.
        ----------------------------------------------------------------------
        costs:
        This is a dictionary containing the cost values per arc.
        ----------------------------------------------------------------------
        time_limit:
        The time limit for the model. It is a positive
        integer, telling the model how much time it maximally has to find a
        solution.
        ----------------------------------------------------------------------
        dummy_nodes:
        A list containing two elements, both being booleans. The first
        element is True if there is a dummy source node in the randomly
        generated graph, the second element is True if there is a dummy sink
        node. This attribute is especially important when generating random
        graphs.
        ----------------------------------------------------------------------
        initial_cost:
        Each time a constraint is added (apart from the time-
        constraint), a new fake-variable is added with a high cost. For this
        reason, it is important to set a high initial cost. The default value
        is 1000.
        ----------------------------------------------------------------------
        objective:
        This is a string, containing either "MAXIMIZE" or "MINIMIZE".
        It tells the subproblem what to do with the modified costs. The
        default value is "MINIMIZE".
        """

        self.arcs = arcs
        self.times = times
        self.costs = costs
        self.dummies = dummy_nodes
        self.time_limit = time_limit
        self.initial_cost = initial_cost
        self.objective = objective

        # Initialising the RMP model:
        self.model = Model("RMP")
        # Initiating the submodel outside its method, because it is used in a
        # loop in the column generation:
        self.sub_model = Model("Subproblem")

        # Create the neighbours for each arc (to be used in the submodel):
        in_neighbours = {}
        out_neighbours = {}
        for i, j in self.arcs:
            if i not in out_neighbours:
                out_neighbours[i] = []
            out_neighbours[i].append(j)
            if j not in in_neighbours:
                in_neighbours[j] = []
            in_neighbours[j].append(i)

        # Initialise the flow variables
        self.x = self.sub_model.addVars(self.arcs, vtype=GRB.CONTINUOUS, ub=1,
                                        lb=0, name="x")

        # Add flow constraints for source, sink and intermediary nodes:
        # Starting with the source and sink nodes. If the graph is randomly
        # generated, the first and last nodes will not necessarily be the
        # source and sink. Thus, if self.dummies is not [False, False], the
        # indices of the source and sink still have to be found. In other
        # words, the source needs to be found by checking if a node is not in
        # the in_neighbours dict. The opposite is used for the sink node.

        for key in out_neighbours.keys():
            if key not in in_neighbours:
                source_idx = key
                self.sub_model.addConstr(quicksum(self.x[source_idx, j]
                                                  for j in
                                                  out_neighbours[source_idx])
                                         == 1, name="Source_Node")
        for key in in_neighbours.keys():
            if key not in out_neighbours:
                sink_idx = key
                self.sub_model.addConstr(quicksum(self.x[i, sink_idx] for i in
                                                  in_neighbours[sink_idx])
                                         == 1, name="Sink_Node")
        print("Source node:", source_idx)
        print("Sink node:", sink_idx)
        for i in range(min(in_neighbours.keys()), len(in_neighbours)):
            if i != sink_idx and i != source_idx:
                self.sub_model.addConstr(quicksum(self.x[i, j] for j in
                                                  out_neighbours[i]) -
                                         quicksum(self.x[j, i] for j in
                                                  in_neighbours[i]) == 0,
                                         name=f"Node_{i}")

    # Next step in this class: checking all attributes for the correct input,
    # using decorators.
    @property
    def arcs(self) -> list:
        """ Returns the arcs of the network. """
        return self._arcs

    @arcs.setter
    def arcs(self, arcs: list) -> None:
        """ Sets the arcs of the network. """
        if isinstance(arcs, list):
            for arc in arcs:
                if not isinstance(arc, tuple):
                    raise ValueError("Arcs must be tuples.")
                if len(arc) != 2:
                    raise ValueError("Arcs must be tuples of length 2.")
                if not isinstance(arc[0], int) or not isinstance(arc[1], int):
                    raise ValueError("Nodes must be integers.")
                if arc[0] == arc[1]:
                    raise ValueError("Nodes of an arc must be different.")
                if arc[0] < 0 or arc[1] < 0:
                    raise ValueError("Nodes of an arc must be positive "
                                     "integers.")
        else:
            raise ValueError("Arcs must be a list of tuples.")
        self._arcs = arcs

    @property
    def costs(self) -> dict:
        """ Returns the costs of the arcs. """
        return self._costs

    @costs.setter
    def costs(self, values: dict) -> None:
        """ Sets the costs of the arcs. """
        if isinstance(values, dict):
            for key in values.keys():
                if not isinstance(key, tuple):
                    raise ValueError("Keys of the costs dictionary must be "
                                     "tuples.")
                if len(key) != 2:
                    raise ValueError("Keys of the costs dictionary must be "
                                     "tuples of length 2.")
                if key not in self.arcs:
                    raise ValueError(f"Key {key} not found in arcs.")
            for val in values.values():
                if not isinstance(val, (int, float)):
                    raise ValueError("Values of the costs dictionary must be "
                                     "integers or floats.")
        else:
            raise ValueError("Costs must be a dictionary.")
        self._costs = values

    @property
    def times(self) -> dict:
        """ Returns the time values per arc. """
        return self._times

    @times.setter
    def times(self, values: dict) -> None:
        """ Sets the time values per arc. """
        if isinstance(values, dict):
            for key in values.keys():
                if not isinstance(key, tuple):
                    raise ValueError("Keys of the times dictionary must be "
                                     "tuples.")
                if len(key) != 2:
                    raise ValueError("Keys of the times dictionary must be "
                                     "tuples of length 2.")
                if key not in self.arcs:
                    raise ValueError(f"Key {key} not found in arcs.")
            for val in values.values():
                if not isinstance(val, (int, float)):
                    raise ValueError("Values of the times dictionary must be "
                                     "integers or floats.")
        else:
            raise ValueError("Times must be a dictionary.")
        self._times = values

    @property
    def time_limit(self) -> int:
        """ Returns the time limit for the model. """
        return self._time_limit

    @time_limit.setter
    def time_limit(self, value: int) -> None:
        """ Sets the time limit for the model. """
        if not isinstance(value, int):
            raise ValueError("Time limit must be an integer.")
        if value <= 0:
            raise ValueError("Time limit must be a positive integer.")
        self._time_limit = value

    @property
    def dummies(self) -> list:
        """ Returns the dummy nodes. """
        return self._dummies

    @dummies.setter
    def dummies(self, values: list) -> None:
        """ Sets the dummy nodes. """
        if not isinstance(values, list):
            raise ValueError("Dummy nodes must be a list.")
        if len(values) != 2:
            raise ValueError("Dummy nodes must be a list of length 2.")
        for val in values:
            if not isinstance(val, bool):
                raise ValueError("Dummy nodes must be booleans.")
        self._dummies = values

    @property
    def initial_cost(self) -> int:
        """ Returns the initial cost for the model. """
        return self._initial_cost

    @initial_cost.setter
    def initial_cost(self, value: int) -> None:
        """ Sets the initial cost for the model. """
        if not isinstance(value, int):
            raise ValueError("Initial cost must be an integer.")
        if value <= 0:
            raise ValueError("Initial cost must be a positive integer.")
        self._initial_cost = value

    @property
    def objective(self) -> list:
        """ Returns the objective of the model. """
        return self._objective

    @objective.setter
    def objective(self, objective: str = "MINIMIZE") -> str:
        """ Sets the objective of the model. """
        if isinstance(objective, str):
            if objective.upper() != "MINIMIZE" and (objective.upper() !=
                                                    "MAXIMIZE"):
                raise ValueError("Objective must be 'MINIMIZE' or "
                                 "'MAXIMIZE'")
        else:
            raise ValueError("Objective must be either 'MAXIMIZE' or "
                             "'MINIMIZE.")
        self._objective = objective.upper()

    # Checking the input of the Simplex method outside its method to avoid
    # cluttering. This method is called only in the Simplex method, because
    # the column generation uses the Simplex method to start its process.
    def input_checker(self, branches: list = False, operators: list = False,
                      branching_vals: list = False) -> None:
        """ A method that checks the input of the Simplex method (this method
        is called, to make the Simplex less cluttered).
        ======================================================================
        ATTRIBUTES:
        branches:
        An optional variable constaining the branching arcs in list. If the
        value is False, column generation without branching is implemented.
        Each arc is a tuple of length 2, containing the start and end node
        of the arc (positive integers). Even if there is only one branch,
        this should still be inputted as a list. The arcs in the branch should
        be in the arcs list. If branches are inputted, you are required to
        give operators and branching_vals as well.
        To branch on multiple arcs at once, you can make the list have double
        brackets. For example, branching on arc (1, 2) and (1, 3) would be
        inputted as: [[(1, 2), (1, 3)]].
        ----------------------------------------------------------------------
        operators:
        A list containing the operators you wish to use in the branching. An
        operator should be inputted as a string in a list. The following three
        operators are allowed: ["=="], [">="], ["<="].
        ----------------------------------------------------------------------
        branching_vals:
        A list contraining the values of the branching constraints. These
        values should be positive integers. If the value is 1, the
        arc is enforced to be in the path. Use value 0 if you don't want the
        arc in the path. Value 1 cannot be paired with "<=" and 0 with ">=".

        An example-input could look like this:
        branches = [[(1, 3), (3, 2)], [(1, 2)]]
        operators = [">=", "=="]
        branching_vals = [1, 1]
        In this case, the branch is enforcing both arcs (1, 3) and (3, 2) in
        1 constraint, and arc (1, 2) is also enforced in an additional
        constraint.
        """

        if branches:
            if not isinstance(branches, list):
                raise ValueError("Branches must be a list.")
            for branch in branches:
                if not isinstance(branch, (list, tuple)):
                    raise ValueError("Branches must be tuples or a list of "
                                     "tuples.")
                if isinstance(branch, tuple):
                    if len(branch) != 2:
                        raise ValueError("Arcs in branches must be tuples of "
                                         "length 2.")
                    if branch[0] < 0 or branch[1] < 0:
                        raise ValueError("Nodes must be positive integers.")
                if isinstance(branch, list):
                    for bran in branch:
                        if not isinstance(bran, tuple):
                            raise ValueError("Branches must be tuples.")
                        if len(bran) != 2:
                            raise ValueError("Arcs in branches must be tuples "
                                             "of length 2.")
                        if bran not in self.arcs:
                            raise ValueError(f"Arc {bran} not found in arcs.")
                        if bran[0] < 0 or bran[1] < 0:
                            raise ValueError("Nodes must be positive "
                                             "integers.")

            if not operators:
                raise ValueError("Operators must be provided.")
            if not branching_vals:
                raise ValueError("Branching values must be provided.")

            if not isinstance(operators, list):
                raise ValueError("Operators must be a list.")
            if len(operators) != len(branches):
                raise ValueError(
                    "Operators must be the same length as branches.")
            for operator in operators:
                if operator not in ["==", ">=", "<="]:
                    raise ValueError("Operators must be ['=='], ['>='], or "
                                     "['<='].")

            if not isinstance(branching_vals, list):
                raise ValueError("Branching values must be a list.")
            if len(branching_vals) != len(branches):
                raise ValueError(
                    "Branching values must be the same length as branches.")
            for i, val in enumerate(branching_vals):
                if not isinstance(val, (int, float)):
                    raise ValueError("Branching values must be integers or "
                                     "floats.")
                if val < 0:
                    raise ValueError("Branching values must be greater than "
                                     "or equal to 0.")
        else:
            if operators or branching_vals:
                raise ValueError("If branches are not provided, operators and"
                                 " branching values cannot be provided.")

    def simplex(self, branches: list = False, operators: list = False,
                branching_vals: list = False) -> None:
        """ Solve the model using the Simplex method.
        ======================================================================
        ATTRIBUTES:
        branches:
        An optional variable constaining the branching arcs in list. If the
        value is False, column generation without branching is implemented.
        Each arc is a tuple of length 2, containing the start and end node
        of the arc (positive integers). Even if there is only one branch,
        this should still be inputted as a list. The arcs in the branch should
        be in the arcs list. If branches are inputted, you are required to
        give operators and branching_vals as well.
        To branch on multiple arcs at once, you can make the list have double
        brackets. For example, branching on arc (1, 2) and (1, 3) would be
        inputted as: [[(1, 2), (1, 3)]].
        ----------------------------------------------------------------------
        operators:
        A list containing the operators you wish to use in the branching. An
        operator should be inputted as a string in a list. The following three
        operators are allowed: ["=="], [">="], ["<="].
        ----------------------------------------------------------------------
        branching_vals:
        A list contraining the values of the branching constraints. These
        values should be integers; either 0 or 1. If the value is 1, the
        arc is enforced to be in the path. Use value 0 if you don't want the
        arc in the path. Value 1 cannot be paired with "<=" and 0 with ">=".

        An example-input could look like this:
        branches = [[(1, 3), (3, 2)], [(1, 2)]]
        operators = [">=", "=="]
        branching_vals = [1, 1]
        In this case, the branch is enforcing both arcs (1, 3) and (3, 2) in
        1 constraint, and arc (1, 2) is also enforced in an additional
        constraint.
        """
        # Checking the input before proceeding:
        self.input_checker(branches, operators, branching_vals)

        # First, a fake variable is created with a large cost and 0 path time,
        # so a feasible (initial) solution is guaranteed (if the model is
        # feasible).
        # Next, the convexity and time constraints are added, followed by the
        # optional branching constraints.
        # The solution to this initial problem will be the kick-starter for
        # the column generation.
        y0 = self.model.addVar(obj=self.initial_cost, name="y0")
        convexity_constraint = self.model.addConstr(y0 == 1)
        time_constraint = self.model.addConstr(y0 * 0 <= self.time_limit,
                                               name="time_constraint")

        # The branch constraints are needed in the column generation, so they
        # are added to the output of this method.
        # Aditionally, the branching variables are stored for the column
        # column generation.
        branch_constraints = []

        # Adding the optional branching constraints to the model is a bit
        # tricky to implement; I chose to have a list as input,
        # containing the equality or inequality operators as strings. This
        # method checks the operator and adds the corresponding constraint
        # to the model, using the right function.
        def equality(constraint, branching_val):
            branch_var = self.model.addVar(obj=self.initial_cost,
                                           name=f"y{constraint + 1}")

            branch_constraint = self.model.addConstr(branch_var ==
                                                     branching_val,
                                                     name="branch_"
                                                          "constraint_"
                                                          f"{constraint + 1}")
            branch_constraints.append(branch_constraint)

        def greater_than(constraint, branching_val):
            branch_var = self.model.addVar(obj=self.initial_cost,
                                           name=f"y{constraint + 1}")

            branch_constraint = self.model.addConstr(branch_var >=
                                                     branching_val,
                                                     name=f"branch_"
                                                          "constraint_"
                                                          f"{constraint + 1}")
            branch_constraints.append(branch_constraint)

        def lesser_than(constraint, branching_val):
            branch_var = self.model.addVar(obj=self.initial_cost,
                                           name=f"y{constraint + 1}")

            branch_constraint = self.model.addConstr(branch_var <=
                                                     branching_val,
                                                     name=f"branch_"
                                                          "constraint_"
                                                          f"{constraint + 1}")
            branch_constraints.append(branch_constraint)

        # Deciding which constraint to add to the model, based on the input
        # of the user:
        if branches:
            for constraint in range(len(branches)):
                if operators[constraint] == "==":
                    equality(constraint, branching_vals[constraint])
                elif operators[constraint] == "<=":
                    lesser_than(constraint, branching_vals[constraint])
                elif operators[constraint] == ">=":
                    greater_than(constraint, branching_vals[constraint])
                else:
                    raise ValueError("Unknown operator: "
                                     f"{operators[constraint]}\n"
                                     "Please only use '==', '>=' and '<='.")
        # At last, the model is updated after adding all constraints. All
        # values are set to 0, apart from the convexity constraint.
        self.model.update()
        values = [1] + [0] * (len(branches) + 1 if branches else 1)
        constraints = ([convexity_constraint, time_constraint] +
                       branch_constraints)

        return values, constraints

    def sub_problem(self, pi_values: list = [0, 0], branches: list = False):
        """ The subprobem returns a new best path based on the dual variables
        of the RMP solution. It is iteratively called by the column
        generation algorithm.
        ======================================================================
        ATTRIBUTES:
        pi_values:
        A list of at least length 2. The list contains the dual variables'
        values, per constraint in the primal problem.
        ----------------------------------------------------------------------
        branches:
        An optional variable constaining the branching arcs in list. If the
        value is False, column generation without branching is implemented.
        Each arc is a tuple of length 2, containing the start and end node
        of the arc (positive integers). Even if there is only one branch,
        this should still be inputted as a list. The arcs in the branch should
        be in the arcs list. If branches are inputted, you are required to
        give operators and branching_vals as well.
        To branch on multiple arcs at once, you can make the list have double
        brackets. For example, branching on arc (1, 2) and (1, 3) would be
        inputted as: [[(1, 2), (1, 3)]].
        """
        print("SUBPROBLEM")
        print("pi values:", pi_values)

        # Every arc gets a modified cost, based on the constraints added to
        # the model. If there is a branch, check if it is saved as a tuple (1
        # arc) or as a list (multiple arcs). Although I ask the user to always
        # input lists, I still give them the freedom to only use a tuple (if
        # there is only one branching constraint). Then, check if the arc is
        # in the branch and adjust the modified cost.
        # The modified cost is based on the dual variables; if there are 2
        # branching constraints and the current arc in the iteration is in
        # the first constraint, the modified cost is adjusted with the dual
        # variable of the first constraint (the dual variable of the second
        # constraint will be 0 for that arc).
        modified_costs = {}
        for arc in self.arcs:
            if branches:
                for i in range(len(branches)):
                    if isinstance(branches[i], list):
                        if arc in branches[i]:
                            modified_costs[arc] = self.costs[arc] - \
                                                  pi_values[1] * \
                                                  self.times[arc] - \
                                                  pi_values[i + 2]
                            break
                        else:
                            modified_costs[arc] = self.costs[arc] - \
                                                  pi_values[1] * \
                                                  self.times[arc]
                    elif isinstance(branches[i], tuple):
                        if arc == branches[i]:
                            modified_costs[arc] = self.costs[arc] - \
                                                  pi_values[1] * \
                                                  self.times[arc] - \
                                                  pi_values[i + 2]
                            break
                        else:
                            modified_costs[arc] = self.costs[arc] - \
                                                  pi_values[1] * \
                                                  self.times[arc]
            else:
                modified_costs[arc] = self.costs[arc] - pi_values[1] * \
                                      self.times[arc]

        # Objective: minimize reduced cost for the path.
        if self.objective == "MINIMIZE":
            self.sub_model.setObjective(quicksum(modified_costs[arc] *
                                                 self.x[arc] for arc in
                                                 self.arcs) - pi_values[0],
                                        GRB.MINIMIZE)
        elif self.objective == "MAXIMIZE":
            self.sub_model.setObjective(quicksum(modified_costs[arc] *
                                                 self.x[arc] for arc in
                                                 self.arcs) - pi_values[0],
                                        GRB.MAXIMIZE)

        # Solve the subproblem:
        self.sub_model.optimize()

        # Extract the path if the reduced cost is negative (for minimization)
        # or positive (for maximization).
        if abs(self.sub_model.ObjVal) < 1.e-6:
            return None, None
        else:
            new_path = [arc for arc in self.arcs if self.x[arc].X > 0.5]
            return new_path

    def column_generation(self, branches: list = False,
                          operators: list = False,
                          branching_vals: list = False) -> None:
        """ Implements the column generation and returns the path(s) found.
        ======================================================================
        ATTRIBUTES:
        branches:
        An optional variable constaining the branching arcs in list. If the
        value is False, column generation without branching is implemented.
        Each arc is a tuple of length 2, containing the start and end node
        of the arc (positive integers). Even if there is only one branch,
        this should still be inputted as a list. The arcs in the branch should
        be in the arcs list. If branches are inputted, you are required to
        give operators and branching_vals as well.
        To branch on multiple arcs at once, you can make the list have double
        brackets. For example, branching on arc (1, 2) and (1, 3) would be
        inputted as: [[(1, 2), (1, 3)]].
        ----------------------------------------------------------------------
        operators:
        A list containing the operators you wish to use in the branching. An
        operator should be inputted as a string in a list. The following three
        operators are allowed: ["=="], [">="], ["<="].
        ----------------------------------------------------------------------
        branching_vals:
        A list contraining the values of the branching constraints. These
        values should be integers; either 0 or 1. If the value is 1, the
        arc is enforced to be in the path. Use value 0 if you don't want the
        arc in the path. Value 1 cannot be paired with "<=" and 0 with ">=".

        An example-input could look like this:
        branches = [[(1, 3), (3, 2)], [(1, 2)]]
        operators = [">=", "=="]
        branching_vals = [1, 1]
        In this case, the branch is enforcing both arcs (1, 3) and (3, 2) in
        1 constraint, and arc (1, 2) is also enforced in an additional
        constraint.
        """

        # Get a first solution, so the column generation can start. No need to
        # check for input values, as the simplex method already does that.
        values, constraints = self.simplex(branches, operators,
                                           branching_vals)
        lambda_vars = {}  # A dictionary that will contain every path found.
        path_times = {}  # A dictionary that will contain the total time of
        # every path found.

        iterations = 0
        while iterations < 30:
            # Solve the RMP:
            self.model.optimize()

            if self.model.status != GRB.OPTIMAL:
                break

            # Get dual variables from RMP:
            pi_values = []
            for index, constraint in enumerate(constraints):
                exponent = 0
                if index >= 2:
                    if branching_vals[index - 2] == 0:
                        exponent = 1
                pi_values.append(constraint.Pi * (-1) ** exponent)
            # If a constraint's value is 0 (meaning it aims to omit), the
            # corresponding dual variable's sign is switched. The first two
            # constraints are always respectively the convexity constraint
            # and the time constraint.

            # Solve subproblem to get a new path and its reduced cost:
            new_path = self.sub_problem(pi_values, branches)

            # The subproblem returns [None, None] as a path if the reduced 
            # cost is non-negative for minimization or non-positive for 
            # maximization problems.
            if new_path[0] is None:
                break

            # Add new path to the RMP with corresponding lambda variable:
            path_cost = sum(self.costs[arc] for arc in new_path)
            path_time = sum(self.times[arc] for arc in new_path)
            values[1] = path_time
            # The second constraint is the time constraint.
            for val in range(2, len(values)):
                values[val] = 0

            # Every time an arc, which is being branched on, is in the path,
            # the constraints value changes (multiple arcs can make this value
            # higher than 1).
            if branches:
                for arc in new_path:
                    for branch in range(len(branches)):
                        if isinstance(branches[branch], list):
                            if arc in branches[branch]:
                                values[branch + 2] += 1
                        elif isinstance(branches[branch], tuple):
                            if arc == branches[branch]:
                                values[branch + 2] += 1

            column = Column(values, constraints)
            lambda_var = self.model.addVar(obj=path_cost, column=column,
                                           vtype=GRB.CONTINUOUS,
                                           name=f"lambda_{len(lambda_vars)}")

            self.model.update()

            if self.dummies[0]:
                nodes = set()
                for arc in self.arcs:
                    nodes.add(arc[0])
                if new_path[0][0] != max(nodes):
                    # If the first arc in the new path is not the dummy
                    # source node, the dummy source node is moved to the
                    # first arc in the new path. This is only done because it
                    # is more readable for the user (the dummy nodes are 
                    # always added as the last node indices in the graph).
                    for i in range(1, len(new_path)):
                        if new_path[i][0] == max(nodes):
                            temp = new_path[i]
                            new_path.pop(i)
                            new_path.insert(0, temp)
                            break

            lambda_vars[lambda_var] = new_path
            path_times[lambda_var] = path_time
            iterations += 1
            if iterations == 30:
                return 0

        # Printing and returning the solution:
        if self.model.status != GRB.INFEASIBLE:
            # Initializing a dictionary with arc flows put to 0 (for
            # displaying individual arc nodes):
            arc_flows = {arc: 0.0 for arc in self.arcs}

            print("\nFinal path's arc flows:")
            for lambda_var, path in lambda_vars.items():
                # Retrieve the variable from the model:
                model_var = self.model.getVarByName(lambda_var.VarName)
                if model_var is not None and model_var.X > 1e-6:
                    path_flow = model_var.X
                    print(f"Path {model_var.VarName} (flow = {path_flow}): "
                          f"{path}")
                    # Add this path's flow to each arc it contains:
                    for arc in path:
                        if arc in arc_flows:
                            arc_flows[arc] += path_flow

            print("\nFinal Individual Arc Flows (non-zero):")
            non_zero_arcflows = {}
            for arc, flow in arc_flows.items():
                if flow > 1e-6:
                    print(f"Arc {arc}: Flow = {flow}")
                    non_zero_arcflows[arc] = flow

            return lambda_vars


##############################################################################
# The times and costs from the Primer:
def primer():
    """ Return the costs and times stated in the Primer in a dict. """
    costs = {(1, 2): 1, (1, 3): 10, (3, 2): 1, (2, 4): 1, (2, 5): 2,
             (3, 5): 12, (3, 4): 5, (4, 5): 10, (4, 6): 1, (5, 6): 2}
    times = {(1, 2): 10, (1, 3): 3, (3, 2): 2, (2, 4): 1, (2, 5): 3,
             (3, 5): 3, (3, 4): 7, (4, 5): 1, (4, 6): 7, (5, 6): 2}
    return costs, times


##############################################################################
# The different examples from the Primer:
def BB0(arcs, costs, times):
    """ Runs the first branch-and-bound iteration from the Primer. """

    print("-" * 100)
    print("Running the first branch-and-bound iteration from the Primer.")
    print("-" * 100)
    column = ColumnGeneration(arcs, times, costs, initial_cost=100)
    column.column_generation()


def BB1(arcs, costs, times):
    """ Runs the second branch-and-bound iteration from the Primer,
    using x_(13) + x_(32) = 0 as a branch. """

    print("-" * 100)
    print("Running the second branch-and-bound iteration from the Primer, "
          "using x_(13) + x_(32) = 0 as a branch")
    print("-" * 100)
    column = ColumnGeneration(arcs, times, costs)
    column.column_generation([[(1, 3), (3, 2)]], ["=="], [0])


def BB2(arcs, costs, times, ):
    """ Runs the third branch-and-bound iteration from the Primer,
    using x_(13) + x_(32) >= 1 as a branch. """

    print("-" * 100)
    print("Running the third branch-and-bound iteration from the Primer, "
          "using x_(13) + x_(32) >= 1 as a branch")
    print("-" * 100)
    column = ColumnGeneration(arcs, times, costs)
    column.column_generation([[(1, 3), (3, 2)]], [">="], [1])


def BB3(arcs, costs, times):
    """ Runs the fourth branch-and-bound iteration from the Primer,
    using x_(13) + x_(32) >= 1 and x_(12) = 0 as a branches. """

    print("-" * 100)
    print("Running the second branch-and-bound iteration from the Primer, "
          "using x_(13) + x_(32) >= 1 and x_(12) = 0 as branches")
    print("-" * 100)
    column = ColumnGeneration(arcs, times, costs)
    column.column_generation([[(1, 3), (3, 2)], [(1, 2)]], [">=", "=="],
                             [1, 0])


def BB4(arcs, costs, times):
    """ Runs the second branch-and-bound iteration from the Primer,
    using x_(13) + x_(32) >= 1 and x_(12) >= 1 as branches. Note that the 
    Gurobi solver will prefer to state that a sum of multiple paths is still 
    feasible, instead of stating that the constraints are infeasible. 
    """

    print("-" * 100)
    print("Running the second branch-and-bound iteration from the Primer, "
          "using x_(13) + x_(32) >= 1 and x_(12) = 1 as branches")
    print("-" * 100)
    column = ColumnGeneration(arcs, times, costs)
    column.column_generation([[(1, 3), (3, 2)], [(1, 2)]], [">=", "=="],
                             [2, 1])
    print("Note that the Gurobi solver will prefer to state that a sum of "
          "multiple paths is still feasible, instead of stating that the "
          "constraints are infeasible. This is a built-in mechanism of "
          "Gurobi.")


##############################################################################
# Randomly generated graphs:
def random_graph(number_nodes: int = 15, number_arcs: int = 30,
                 plot: bool = False):
    """ Generates a random directed acyclic graph (DAC). As the rest of this
    code relies on 1 source node and 1 sink node, this function checks if
    there are multiple sources or sinks generated. If so, a dummy-node is
    created with a cost and time-cost of 0. This way, there's always 1
    source node and 1 sink node (and the requested amount of actual nodes
    remains the same).
    ======================================================================
    VARIABLES:
    number_nodes:
    The amount of nodes (as integers) you want to have in the generated
    DAC.
    ----------------------------------------------------------------------
    number_arcs:
    The number of arcs (as integers). This is typically double the amount
    of nodes.
    ----------------------------------------------------------------------
    plot:
    A boolean variable, telling the function if you want the newly
    generated graph plotted or not.
   """
    # First, a random seed is set, for reproducibility:
    random.seed(42)

    # Generating a random undirected graph with a fixed number of edges,
    # without loops (undirected for now, so the mode "acyclic" can be used):
    g = ig.Graph.Erdos_Renyi(n=number_nodes, m=number_arcs, directed=False,
                             loops=False)

    # Converting to a Directed Acyclic Graph (DAG)
    g.to_directed(mode="acyclic")

    # Sometimes, the randomly generated graph contains nodes with no incoming
    # or outgoing arcs. These nodes needs to be deleted (this will
    # unfortunately likely cause the amount of nodes to be less than
    # requested).
    deletes = 0
    for node in g.vs:
        if node.indegree() == 0 and node.outdegree() == 0:
            node.delete()
            deletes += 1
    if deletes > 0:
        print("Deleted", deletes, "nodes with no edges.")

    # Adding dummy nodes and sinks, because a randomly generated graph can
    # have multiple sources and sinks. These nodes are added to make sure
    # there is only one source and one sink in the graph.
    sources = [v.index for v in g.vs if v.indegree() == 0]
    sinks = [v.index for v in g.vs if v.outdegree() == 0]
    dummies = []
    dummy_source_index = 0
    dummy_sink_index = 0

    if len(sources) > 1:
        g.add_vertex("dummy_source")
        # Because a vertex was added, the vcount is now 1 more than before,
        # which means that the dummy source will be the (next to) last node in
        # the graph.
        dummy_source_index = g.vcount() - 1  
        # iGraph starts its nodes at index 0.
        dummies.append(True)
        for source_node in sources:
            g.add_edge(dummy_source_index, source_node)
    else:
        dummies.append(False)
    if len(sinks) > 1:
        g.add_vertex("dummy_sink")
        # The vcount can be reused, because the dummy source is added before
        # the dummy sink.
        dummy_sink_index = g.vcount() - 1
        dummies.append(True)
        for sink_node in sinks:
            if sink_node != dummy_source_index:
                g.add_edge(sink_node, dummy_sink_index)
    else:
        dummies.append(False)
    # Now, because the dummies list contains information on if a dummy node
    # and/or a dummy sink is created, the column generation class handles the
    # source correctly. The dummy source has either the last or next-to-last
    # index, depending on if there is a dummy sink.

    arcs = []
    costs = {}
    times = {}
    dummy_indices = set()
    if dummy_source_index != 0:
        dummy_indices.add(dummy_source_index)
    if dummy_sink_index != 0:
        dummy_indices.add(dummy_sink_index)

    # Assigning random costs and times to arcs that are not from or to the
    # dummy nodes:
    for arc in g.es:
        source = arc.source
        target = arc.target
        arcs.append((source, target))
        if source in dummy_indices or target in dummy_indices:
            cost = 0
            time = 0
        else:
            cost = random.randint(1, 20)
            time = random.randint(1, 20)
        costs[(source, target)] = cost
        times[(source, target)] = time

    print("Graph with", number_nodes, "nodes and", number_arcs,
          "edges generated.")
    print("Costs:", costs, "\n")
    print("Times:", times)
    print("Dummy nodes", dummies)
    if plot:
        edge_labels = [f"({costs[edge.tuple]}, {times[edge.tuple]})"
                       for edge in g.es]
        fig, ax = plt.subplots()
        ig.plot(
            g,
            target=ax,
            layout="sugiyama",
            vertex_size=30,
            vertex_color="grey",
            vertex_label=g.vs.indices,
            vertex_label_size=9,
            vertex_label_dist=0,  # Placing label inside the node.
            vertex_label_color="white",
            edge_color="#222",
            edge_width=1,
            edge_label=edge_labels,
            edge_label_size=5,
            edge_align_label=True,
            edge_background="white",
            title="Randomly generated graph (with added dummy nodes if "
                  "there were multiple sources and/or sinks)",
        )
        plt.show()

    return costs, times, dummies


##############################################################################
# Experimental:
def OneArc(arcs, costs, times):
    """ Runs and times every iteration from the Primer, but with only one arc
    per branch. """
    # BB1:
    start_time = datetime.now()
    print("-" * 100)
    print("Running BB1 with 1 constraint")
    print("-" * 100)
    column = ColumnGeneration(arcs, times, costs)
    column.column_generation([[(1, 3), (3, 2)]], ["=="], [0])
    end_time = datetime.now()

    start_time2 = datetime.now()
    print("-" * 100)
    print("Running BB1 with 2 constraints")
    print("-" * 100)
    column = ColumnGeneration(arcs, times, costs)
    column.column_generation([[(1, 3)], [(3, 2)]], ["==", "=="], [0, 0])
    end_time2 = datetime.now()

    # BB2:
    start_time3 = datetime.now()
    print("-" * 100)
    print("Running BB2 with 1 constraint")
    print("-" * 100)
    column = ColumnGeneration(arcs, times, costs)
    column.column_generation([[(1, 3), (3, 2)]], [">="], [1])
    end_time3 = datetime.now()

    start_time4 = datetime.now()
    print("-" * 100)
    print("Running BB2 with 2 constraints")
    print("-" * 100)
    column = ColumnGeneration(arcs, times, costs)
    column.column_generation([[(1, 3)], [(3, 2)]], [">=", ">="], [1, 1])
    end_time4 = datetime.now()

    # BB3:
    start_time5 = datetime.now()
    print("-" * 100)
    print("Running BB3 with 2 constraints")
    print("-" * 100)
    column = ColumnGeneration(arcs, times, costs)
    column.column_generation([[(1, 3), (3, 2)], [(1, 2)]], [">=", "=="],
                             [1, 0])
    end_time5 = datetime.now()

    start_time6 = datetime.now()
    print("-" * 100)
    print("Running BB3 with 3 constraints")
    print("-" * 100)
    column = ColumnGeneration(arcs, times, costs)
    column.column_generation([[(1, 3)], [(3, 2)], [(1, 2)]], [">=", ">=",
                                                              "=="],
                             [1, 1, 0])
    end_time6 = datetime.now()

    print("=" * 100)
    print("BB1:", end_time - start_time)
    print("BB1 with 1 more constraint", end_time2 - start_time2)
    print("=" * 100)
    print("BB2:", end_time3 - start_time3)
    print("BB2 with 1 more constraint", end_time4 - start_time4)
    print("=" * 100)
    print("BB3:", end_time5 - start_time5)
    print("BB3 with 1 more constraint", end_time6 - start_time6)


def WholePath(arcs, costs, times):
    """ Runs the example from the Primer, but with the whole path as
    branch. """

    print("-" * 100)
    print("Running BB0")
    print("-" * 100)
    # Running the BB0 column generation as the starter to this process:
    column = ColumnGeneration(arcs, times, costs)
    lambda_vars = column.column_generation()
    # Now, the lambda_vars are the variables that are used to branch on: the
    # branch will be the whole path. To do this, every lambda variable that
    # has a fractional value gets looped over. In the loop, there two branches
    # are created: one with the path in the model and one without the path in
    # the model (via constraints). After both newly created models finished, a
    # check is done to see if the model is feasible (if the sum of lambda-
    # values is 1). If not, this is printed alongside the gurobi output.

    # The new variables are stored in a dictionary, with the name of the
    # variable as key and the path as value. This is done to make it easier to
    # print the results of the branching. In the output, the branching
    # paths and values are printed.
    new_vars = {}
    if lambda_vars != 0:
        for var, path in lambda_vars.items():
            if 1 > var.X > 1e-6:
                

                # First branching:
                start_time = datetime.now()
                first_branch = ColumnGeneration(arcs, times, costs)
                print("=" * 100)
                print("Branching (enforcing) on path:", path)
                brancher1 = first_branch.column_generation([path], [">="],
                                                           [1])
                end_time = datetime.now()

                # Collect keys to remove in a separate list to avoid modifying
                # the dictionary during iteration (this gives an error):
                keys_to_remove = [var1 for var1 in brancher1 if var1.X < 1e-6]
                for key in keys_to_remove:
                    brancher1.pop(key)
                new_vars[f"{var.VarName} >=1"] = brancher1

                summ = 0
                for paths in brancher1:
                    summ += paths.X
                if summ < 0.99:
                    print("This model is infeasible.")

                # Second branching:
                start_time2 = datetime.now()
                second_branch = ColumnGeneration(arcs, times, costs)
                print("=" * 100)
                print("Branching (omitting) on path:", path)
                brancher2 = second_branch.column_generation([path], ["=="],
                                                            [0])
                end_time2 = datetime.now()

                keys_to_remove = [var2 for var2 in brancher2 if var2.X < 1e-6]
                for key in keys_to_remove:
                    brancher2.pop(key)
                new_vars[f"{var.VarName} ==0"] = brancher2

                summ = 0
                for paths in brancher2:
                    summ += paths.X
                if summ < 0.99:
                    print("This model is infeasible.")

                print("-" * 100)
                print(f"Branching on {var.VarName} (value {var.X}) - {path} "
                      "results in:\n"
                      "\n>= 1:", new_vars[f"{var.VarName} >=1"], "with a time of:", end_time - start_time,
                      "\n== 0:", new_vars[f"{var.VarName} ==0"], "with a time of:", end_time2 - start_time2)
