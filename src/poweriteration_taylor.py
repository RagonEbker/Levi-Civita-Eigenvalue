#Author Ragon Ebker
#Note: Should only be used for small exampls since its very slow

from sympy import diff, factorial, symbols
import numpy as np

eps = symbols('eps')

def seriesExpansionEps(f,n_terms):
  g = 0
  for i in range(0,n_terms):
    g = g + (diff(f,eps,i).evalf(subs={eps: 0}))/factorial(i)*eps**i
  return g

def normalize_lc_2(v,n_terms):

  dim = v.shape[0]
  norm_b_k = norm_lc_2(v,dim)
  for i in range(dim):
     v[i] =seriesExpansionEps(v[i]/norm_b_k,n_terms)
  return v

def norm_lc_2(v,dim):
   norm = 0
   
   for i in range(0,dim):
      norm = norm + (v[i]**2)
   return seriesExpansionEps(norm**(1/2),5)

def rayleigh_coeff_lc(A,v):
  return v.T @ (A @ v)

def power_iteration(A: np.array, num_iterations: int,start_vec : np.array,n_terms):

  b_k = start_vec
  
  for i in range(num_iterations):
      print("Step " + str(i))
      b_k = A @ b_k

      b_k = normalize_lc_2(b_k,n_terms)

  return b_k