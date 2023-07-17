#Author Ragon Ebker


from sympy import RR
import numpy as np
from sympy.polys.ring_series import ring, rs_nth_root, rs_series_inversion, rs_trunc

R, eps = ring('eps', RR)

def normalize_lc_2(v : np.array,n_s : int):
  dim = v.shape[0]
  norm_b_k = norm_lc_2(v,dim,n_s)
  n_invers = rs_series_inversion(norm_b_k,eps,n_s)
  for i in range(dim):
     v[i] = rs_trunc(v[i] * n_invers,eps,n_s)
  return v

def norm_lc_2(v : np.array,dim : int,n_s : int):
   norm = 0
   
   for i in range(0,dim):
      norm = norm + (v[i]**2)
   return rs_nth_root(norm,2,eps,n_s)

def rayleigh_coeff_lc_wo_norm(A : np.array,v : np.array,n_s):
  return rs_trunc((v.T @ (A @ v)),eps,n_s)

def power_iteration_lc(A: np.array,start_vec : np.array, num_iterations: int, n_s : int):
  """
  Calculates an eigenvector from the eigenspace of the dominant eigenvalue iteratively
  
  Parameters
  ----------
  A             : Levi-Civita Matrix whos eigenvectors we want to calculate
  start_vec     : Start vector for the power iteration, should have a non-vanishing component
                  in the dominant eigenspace (this is numerically negligible)
  num_iterations: Number of iterations
  n_s           : Number of epsilon terms we want conserve in each iteration

  Return:
  -------
  b_k : An approximation of our eigenvector
  """
  b_k = start_vec
  
  for i in range(num_iterations):
      print("Step: " + str(i))
      b_k = A @ b_k

      b_k = normalize_lc_2(b_k,n_s)
      
  return b_k

def power_iteration_lc_w_steps(A: np.array, num_iterations: int,start_vec : np.array, n_s : int):
  
  b_k = start_vec
  b_ks = []
  for i in range(num_iterations):
      print("Step: " + str(i))
      b_k = A @ b_k

      b_k = normalize_lc_2(b_k,n_s)
      b_ks.append(b_k)
  return b_ks