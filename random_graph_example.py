#Author Ragon Ebker

from src.create_random_graph import create_graph_example

#dimension of the matrix
dim = 3
#Max Power of epsilon
n_eps = 2
#max number for the random coefficients
max_n = 10
#number of epsilon terms for representation of lc number
n_s = 10

create_graph_example(dim,n_eps,max_n,n_s)

#Example for creating multiple graph
# max_dim = 31
# max_eps = 10
# create_multiple(max_dim,max_eps,n_s)