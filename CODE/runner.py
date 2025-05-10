# Author: Alexander Colman, student Physics at the University of Antwerp
# Contact: Alexander dot colman1 AT gmail dot com
# Latest update: 10/05/2025

from datetime import datetime
import random

from implementation import ColumnGeneration, primer, \
    BB0, BB1, BB2, BB3, BB4, OneArc, WholePath, random_graph


# This script runs the various methods and functions defined in
# implementation.py

costs, times = primer()
arcs = list(times.keys())

# Run BB0:
# BB0(arcs, costs, times)

# Run BB1:
# BB1(arcs, costs, times)

# Run BB2:
# BB2(arcs, costs, times)

# Run BB3:
# BB3(arcs, costs, times)

# Run BB4:
# BB4(arcs, costs, times)

# Run OneArc:
# OneArc(arcs, costs, times)

# Run WholePath:
# WholePath(arcs, costs, times)

# Run a randomly generated graph:
# num_nodes = 15
# costs, times, dummies = random_graph(num_nodes, num_nodes * 2, True)
# arcs = list(times.keys())
# branc = []
# for _ in range(2):
#     branc.append(arcs[random.randint(0, len(arcs) - 1)])
# print("Branching on:", branc)

# kolom = ColumnGeneration(arcs, times, costs, num_nodes * 2, dummies, 1000, 
#                          "MINIMIZE")
# kolom2 = ColumnGeneration(arcs, times, costs, num_nodes * 2, dummies, 1000, 
#                          "MINIMIZE")

# start_time = datetime.now()
# kolom.column_generation([branch for branch in branc], [">="] * len(branc), 
#                         [1] * len(branc))
# end_time = datetime.now()

# print("Time taken for column generation with all branches in other "
#       "constraints:", end_time - start_time)

# start_time2 = datetime.now()
# kolom.column_generation([branc], [">="], [1])
# end_time2 = datetime.now()

# print("Time taken for column generation with all branches in other "
#       "constraints:", end_time - start_time)
# print("Time taken for column generation with all branches in one list:",
#       end_time2 - start_time2)
