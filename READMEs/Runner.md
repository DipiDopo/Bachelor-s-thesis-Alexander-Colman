runner.py runs the code in implementation.py for various cases.

Initially, almost all code is commented, but uncommenting a line of code will make the 
file run the column generation algorithm, as implemented in implementation.py. 

The first four codes run the example from the Primer. 

The next two run the OneArc and WholePath branching researches (performed on the Primer
model as well). These offer some results, but the length of the paths is too short to 
get significant time differences between the results; the relative time difference is 
always less than 0.05s.

The last code is used to generate a random graph (DAG), with randomly assigned costs and 
time-costs. The computation time research on this generated graph is not that useful
though, as the found paths are typically very short (even in 500 node-networks) and/or
have little branching options. This means the same problem arises as when trying to time 
branching tactics in the Primer model.

Taking some inspiration from the Primer model, a new model was constructed, and its path
length was always 100 arcs. This graph is used in the thesis and can be found in 
research.py.
