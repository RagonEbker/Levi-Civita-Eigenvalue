#Author Ragon Ebker
#Note: Should only be used for small examples because its very slow

from src.poweriteration_taylor import power_iteration,seriesExpansionEps,rayleigh_coeff_lc
import numpy as np
from sympy import symbols

eps = symbols('eps', real = True)

#example matrix
A = np.array([[1,1+eps],[1-eps,1]])

n_iterations = 10
dim = A.shape[0]
start_vec = np.random.rand(dim)
n_terms = 3

result = power_iteration(A,n_iterations,start_vec,n_terms)
ew = seriesExpansionEps(rayleigh_coeff_lc(A,result),n_terms)

print(ew)