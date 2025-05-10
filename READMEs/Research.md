As all previous graphs failed to produce a meaningful outcome for the computation times
of different branching constraints, a new model was created in research.py.

A function create_LP creates an LP with 100 layers of 2 nodes. Every node is connected to
the next two nodes, and the graph contains only 1 sink and 1 source node. Every arc is 
given a random cost and time-cost between 1 and 20. Afterwards, this graph's costs and 
time-costs are saved to a text file, so the solutions are reproducable. The create_LP 
function is optional, as the graph used in the thesis is already save and added as
text file to this GitHub repository.

Next, the costs and time-costs can be read from the text file, after which an initial 
solution is calculated. In optimal solution of the saved model, there are 4 non-integer 
arcs. 

The time_branches function times the calculations of different branching strategies. The
three implemented strategies are: constrainting two arcs in one constraint, dividing
the same two arcs into two separate constraints, and branching a whole path, with every
arc in one constraint. The function generates the same ColumnGeneration object three 
times, because if the same model is re-run with different branches, the previous
constraints are still stored. This would lead to an incorrect (or infeasible) solution.

time_branches is called 50 times, and the computation time of every tactic is stored
per iteration. Each tactic's time is then plotted iteration, alongside the mean and 
error on that mean. To only plot the whole path, comment every plotting line apart from the
one with label "Whole path" and uncomment every line below that code section.
