# Author: Alexander Colman, student Physics at the University of Antwerp
# Contact: Alexander dot colman1 AT gmail dot com
# Latest update: 10/05/2025


from datetime import datetime
import matplotlib.pyplot as plt
import networkx as nx
import random
from scipy.stats import sem

from implementation import ColumnGeneration


# This file researches the computation time of three branching tactics on a 
# self-made graph.


def create_LP():
    # Number of steps between source and sink:
    levels = 100  

    # Initialize graph
    G = nx.DiGraph()

    # Create a sink and a source node:
    source = 0
    sink = 2 * levels + 1
    G.add_node(source)
    G.add_node(sink)

    # Create intermediate nodes (every node is connected to the next two):
    level_nodes = {}
    node_id = 1
    for level in range(levels):
        level_nodes[level] = (node_id, node_id + 1)
        G.add_node(node_id)
        G.add_node(node_id + 1)
        node_id += 2

    # Add edges from source to first two nodes:
    for node in level_nodes[0]:
        G.add_edge(source, node)

    # Connecting all intermediary nodes:
    for level in range(levels - 1):
        current_nodes = level_nodes[level]
        next_nodes = level_nodes[level + 1]
        for u in current_nodes:
            for v in next_nodes:
                G.add_edge(u, v)

    # Add edges from last 2 nodes to sink:
    for node in level_nodes[levels - 1]:
        G.add_edge(node, sink)

    # Plotting:
    plt.figure(figsize=(12, 6))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(G, pos, with_labels=True, node_size=500, arrows=True)
    plt.title("Binary Layered DAG")
    plt.show()

    # Assign random costs and time-costs:
    for u, v in G.edges():
        G[u][v]['cost'] = random.randint(1, 20)
        G[u][v]['time_cost'] = random.randint(1, 20)

    # Assign the random costs and time-costs to the arcs in the graph:
    costs= {}
    times = {}
    for u, v in G.edges():
        costs[(u, v)] = G.edges[u, v]['cost']
        times[(u, v)] = G.edges[u, v]['time_cost']

    # Save costs and times to a text file:
    with open("graph_data.txt", "w") as f:
        f.write("COSTS\n")
        for (u, v), cost in costs.items():
            f.write(f"{u},{v}:{cost}\n")
        f.write("TIMES\n")
        for (u, v), time in times.items():
            f.write(f"{u},{v}:{time}\n")
    print("Costs and times saved to 'graph_data.txt'.")

    # Summarize the graph:
    print(f"Graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} "
          "edges.")
    print(f"Source: {source}, Sink: {sink}")
    print(f"Shortest path length: {nx.shortest_path_length(G, source, sink)}")


# Create an LP with random costs and times (optional; the currently saved LP 
# is the one used in the thesis):
# create_LP()

# Load saved costs and times from file:
costs = {}
times = {}
mode = None

with open("graph_data.txt", "r") as f:
    for line in f:
        line = line.strip()
        if line == "COSTS":
            mode = "costs"
        elif line == "TIMES":
            mode = "times"
        elif line and mode:
            coords, value = line.split(":")
            u, v = map(int, coords.split(","))
            value = int(value)
            if mode == "costs":
                costs[(u, v)] = value
            elif mode == "times":
                times[(u, v)] = value
arcs = list(times.keys())

# Initial solving of the model, to check the solution (which paths are found,
# what their flows are, and which arcs have non-integer flow):
column = ColumnGeneration(arcs, times, costs, 1000, initial_cost=10000)
lambda_vars = column.column_generation()

# Creating a list which contains all non-integer paths:
branches = []
for var, path in lambda_vars.items():
    if 1 > var.X > 0.01:
        branches.append(path)

# Initialzing lists for storing the computation times:
two_constr = []
one_constr = []
whole_path = []

def time_branches(symbol: str = ">=", value: int = 1, value2: int = 2, 
                  value_path: int = 100, path_index: int = 1):
    """ 
    Times the calculations of different branching strategies.
    ==========================================================================
    ATTRIBUTES:
    symbol:
    A string. This can either be ">=" or "==" and is part of determining if 
    the branching will be enforced or omitted.
    --------------------------------------------------------------------------
    value:
    An integer, which can only be 0 or 1. This is the constraint value for the
    2 two-arc branch tactics.
    --------------------------------------------------------------------------
    value2:
    An integer, which can only be 0 or 2. This is the constraint value for the
    2 two-arc branch tactics.
    --------------------------------------------------------------------------
    value_path:
    An integer between 0 and 100. This is the constraint value for branching 
    a whole path.
    --------------------------------------------------------------------------
    path_index:
    An integer, which can only be 0 or 1. This is the index of the path which
    is being branched on. Index 0 refers to the path with the least flow; 
    index 1 refers to the path with the most flow.       
    """

    if symbol == ">=":
       pass 
    elif symbol == "==":
        pass
    else:
        raise ValueError('The branching symbol can only be ">=" or "==".')

    if value == 0: 
       pass
    elif value == 1:
        pass
    else:
        raise ValueError("The constraint value for the two-arcs-two-"
                         "constraints branching can only be between 0 and 1.")

    if value2 == 0: 
       pass
    elif value2 == 2:
        pass
    else:
        raise ValueError("The constraint value for the two-arcs-one-"
                         "constraint branching can only be between 0 and 2.")

    if value_path not in range(101):
        ValueError("The constraint value for branching a whole path can only "
                   "be between 0 and 100.")
    if path_index == 0:
        pass 
    elif path_index == 1:
        pass 
    else:
        ValueError("The index of which path to branch on can only be 0 and "
                   "1.")

    # If the same model is re-run with different branches, the previous
    # constraints are still stored. This would lead to an incorrect (or 
    # infeasible) solution.
    column2 = ColumnGeneration(arcs, times, costs, 1000, initial_cost=10000)
    column3 = ColumnGeneration(arcs, times, costs, 1000, initial_cost=10000)
    column4 = ColumnGeneration(arcs, times, costs, 1000, initial_cost=10000)
    
    # Timing the 3 different branching tactics; the two-arcs branches are 
    # chosen to be index 32 and 33 in the path (this was hand-picked):
    print("=" * 100)
    print(f"branching on the following arcs: {branches[path_index][32:34]} "
          "-> all arcs in seperate constraints")
    start_time1 = datetime.now()
    column2.column_generation(branches[path_index][32:34], [symbol] * 2, 
                              [value] * 2)
    end_time1 = datetime.now()

    print("=" * 100)
    print(f"branching on the following arcs: {branches[path_index][32:34]} "
          "-> all arcs in the same constraint")
    start_time2 = datetime.now()
    column3.column_generation([branches[path_index][32:34]], [symbol], 
                              [value2])
    end_time2 = datetime.now()

    print("=" * 100)
    print("branching on the whole path:", branches)
    start_time3 = datetime.now()
    column4.column_generation([branches[path_index]], [symbol], [value_path])
    end_time3 = datetime.now()

    print("Time taken for column generation with all branches in different "
         "constraints:", end_time1 - start_time1)
    delta = end_time1 - start_time1
    two_constr.append(delta.total_seconds())
    print("Time taken for column generation with all branches in the same "
          "constraint:", end_time2 - start_time2)
    delta = end_time2 - start_time2
    one_constr.append(delta.total_seconds())
    print("Time taken for column generation with the whole path in one "
          "constraint:", end_time3 - start_time3)
    delta = end_time3 - start_time3
    whole_path.append(delta.total_seconds())

# Timing the branching tactics 50 times:
for _ in range(50):
    time_branches(">=", 1, 2, 100, 1)

with open("graph_times.txt", "w") as f:
    f.write("OneConstr, TwoConstr, WholePath\n")
    for time in range(len(one_constr)):
        f.write(f"{one_constr[time]} {two_constr[time]} {whole_path[time]}\n")

# Calculating the mean and standard error of the times:
avg_two = sum(two_constr) / len(two_constr)
err_two = sem(two_constr)
avg_one = sum(one_constr) / len(one_constr)
err_one = sem(one_constr)
avg_whole = sum(whole_path) / len(whole_path)
err_whole = sem(whole_path)
print("Average time taken for column generation with both branches in "
      "different constraints:", avg_two, "with an error of", err_two)
print("Average time taken for column generation with both branches in the "
      "same constraint:", avg_one, "with an error of", err_one)
print("Average time taken for column generation for branching the whole "
      "path:", avg_whole, "with an error of", err_whole)

# Plotting the calculation time, the mean and the error on that mean per 
# tactic:
plt.plot(range(1,51), two_constr, color="b" , label="Two constraints")
plt.plot(range(1,51), one_constr, color="g", label="One constraint")
plt.plot(range(1,51), whole_path, color="r", label="Whole path")
plt.hlines(y=[avg_two, avg_one, avg_whole], xmin=1, xmax=51, 
           colors=['b', 'g', 'r'], linestyles=["--", "--", "--"], 
           label="Averages")
plt.hlines(y=[avg_two + err_two, avg_one + err_one, avg_whole + err_whole], 
           xmin=1, xmax=51, colors=['b', 'g', 'r'], 
           linestyles=[":", ":", ":"], label="Errors")
plt.hlines(y=[avg_two - err_two, avg_one - err_one, avg_whole - err_whole], 
           xmin=1, xmax=51, colors=['b', 'g', 'r'], 
           linestyles=[":", ":" , ":"])

# To only plot the whole path (comment every line apart from the one with 
# label "Whole path"):
# plt.hlines(y=[avg_whole], xmin=1, xmax=50, colors=['r'], linestyles=["--"], 
#            label="Average")
# plt.hlines(y=[avg_whole + err_whole], xmin=1, xmax=50, colors=['r'], 
#            linestyles=[":"], label="Error")
# plt.hlines(y=[avg_whole - err_whole], xmin=1, xmax=50, colors=['r'], 
#            linestyles=[":"])

plt.xlabel("Iteration", fontsize=16)
plt.xlim(1, 50)
plt.xticks(range(1, 51))
plt.ylabel("Duration (s)", fontsize=16)
plt.title("Computation time for different branching tactics", fontsize=16)
plt.legend(loc="lower right", fontsize=16)
plt.show()
