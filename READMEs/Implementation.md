implementation.py contains the code for running the Simplex method, column generation
and Branch-and-Bound algorithms.

This text file gives an overview of what can be found in the Python file.

To start, a class ColumnGeneration is created. The attributes of this class and their
respective description can be found in the docstring of the init method (a further 
explanation on their use is found in the section about the Simplex algorithm and 
using randomly generated graphs). The default values of these attributes are chosen 
based on the most frequent inputs needed to solve the example model found in the 
Primer. The init function already assigns the flow constraints from the source to 
the sink node; this has to be done for the subproblem, but assigning these 
constraints in the subproblem method, will cause the code to malfunction. This is 
because the subproblem method is iteratively called in the column generation process, 
which would make it so new flow constraints are added every iteration of the column 
generation.

The next method in this class is the input_checker. As can be deducted from its name, 
it checks the inputs of the Simplex algorithm. These checks are done in a seperate 
method, as it would otherwise make the simplex method cluttered and hard to read. 
The variables it checks are the branches (arcs, stored in a tuple, which need to be 
inside a list), operators for these branches (e.g. ">=") and the branching values. 
The branches list should be a list containing one or more arcs, depending on the 
branching you would like to perform. If you want multiple arcs inside the same 
branching constraint, you should make a list of this branching list. To be more 
specific: if you want to branch on (1, 3) and (3, 2) in separate constraints, your 
branches variable should look like this: [(1, 3), (3, 2)], but if you want both arcs 
to be used in the same constraint, the branches variable should look like this: 
[[(1, 3), (3, 2)]]. The operators variable should be a list , which can only contain 
the following strings: ">=", "==" and "<=". The amount of constraints given in 
branches should be the same as the length of the operators variable: each constraint 
gets their own operator. Same goes for the branching values: each constraint gets 
their own value. An example input would be: branches = [(1, 3), (3, 2)], 
operators = ["==", "=="], branching_vals = [1, 0]. This would make a constraint on 
the arc (1, 3), which needs to have a flow of 1, while the arc (3, 2) needs to have 
a flow of 0.

The simplex method uses the previously described input_checker. After every input 
is checked and no errors are generated, a fake variable is created and two constraints
are already made, as described in the Primer. Next, the simplex method checks which 
operator was given to a constraint (if there are constraints); a corresponding 
constraint is then added according to which branch what branching is inputted. 
Finally, every constraint and their values are returned, so they can be used in the 
column_generation method.

The subproblem solves creates a new best path, based on the dual variables generated 
in the RMP (which is solved in the column_generation method). Its inputs are, of 
course, the dual variables of the RMP, but also the branches (the same list as in the
simplex method). To calculate the modified costs for each arc, the subproblem method 
has to know which branches are implemented. Based on these branches and their dual
variables, a path is chosen based on the reduced costs of every arc.

At last, the column_generation method implements the column generation and returns
the optimal path(s) found. Its inputs are the same as in the simplex method. Firstly,
the simplex method is called, not only to get the constraints and their values, but 
also to check if the inputs of this method are valid. Next, the column generation
algorithm can commence: the RMP is solved, after which the dual variables get 
extracted (if a constraint is to omit (a) certain arc(s), the sign in front of its 
value should switch. After this, the subproblem is called, using these dual variables,
to generate a new path. Afterwards, a new variable gets created with the values of the
constraints adapted to how much of the arcs in a constraint are in the new path (if 
there are multiple arcs in a constraint). To end the iteration, this new variable gets
stored in a dictionary, along with its total path time.

After the class, a lot of functions were created; these can be called in runner.py, so
that file stays readable and clean.

Firstly, the times and costs can be called which are used in the example model of the
Primer.

Next, every single Branch-and-Bound iteration of the Primer is implemented. These can
be called to get the same solution as in the Primer.

Thirdly, a function was created to randomly generate graphs, with a specified number 
of nodes and arcs. This new graph is a DAG (directed, acyclic graph), so the 
time-constraint holds its meaning. If the generated graph contains multiple sinks 
and/or sources, dummy variables can be created. These are added to the model and 
connect to all original sinks or nodes, so the graph only has one sink and one source. Their costs and time-costs are set to 0, so there was no net change in the model. All other arcs get a random cost and time-cost between 1 and 20. After generating, the user can optionally tell the function to plot the graph.

The OneArc function was created to compare the computation times between different 
branching tactics in the Primer.

The WholePath also uses the Primer's model: the BB0-solution has two non-integer
paths, so this function branches on every arc in one of the paths at once.

Although the previous three functions all tried to determine a relation between the 
computation time and a branching tactic, they all had problems with their path length. The model in the Primer only has a maximum path length of 5 arcs, which means the computation time will always be lower than 0.05s. A similar issue arrises when randomly generating a graph; although its path may be long (which isn't always the case, even when using 500 nodes), there aren't always many possible branching arcs (arcs with non-integer values). Even when importing a large workflow model with 1000 nodes and more than 2000 arcs, the path that was always found was of length 6. Taking some inspiration from the Primer model, a new model was constructed, and its path length was always 100 arcs. This graph is used in the thesis and can be found in research.py.
