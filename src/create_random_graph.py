#Author Ragon Ebker

import numpy as np
import random
import networkx as nx
from sympy import RR
from sympy.polys.ring_series import PolyElement, ring
from sympy.polys.ring_series import _series_inversion1
import os

R, eps = ring('eps', RR)

def create_random_lc_and_real(dim : int, n_eps: int,  max_n: int):

  """
  Create random levi civita matrix A and its real part R_A

  Parameters
  ----------
  dim : int
      dimension of the matrix/graph
  n_eps : int
      maximum power of epsilon
  max_n : int
      maximum number for random coefficients
  
  

  Returns
  -------
  A : ndarray
    LC Matrix (dim,dim) 
  A_R : 
    Real part of the LC Matrix

  """

  #Full Matrix with epsilon parts
  A = np.zeros((dim,dim),PolyElement)

  #Only real part
  R_A = np.zeros((dim,dim),float)

  #prob: probablity that an eps^i, i in [0,eps] arrives in a term. 
  r1 = 1

  for i in range(0,dim):
    for j in range(0,dim):
      if (j > i):
        r = random.randint(0,1)
        if(r > 0):  
          
          a = 0 
          for k in range(0,n_eps+1):
              d = random.randint(r1,1)
              if (d> 0):
                b = random.randint(0,max_n)

                a += b* eps**(k)
                if (k==0):
                  R_A[i,j] = (float) (b)
                  R_A[j,i] = (float) (b)
          A[i,j] =  a
          A[j,i] =  a
  return A,R_A

def laplacian(A : np.array,is_lc,n_s):
  #Creates the laplacian of matrix

  #A: input matrix

  #returns the laplacian of A
  #dimension of our Matrix
  n = A.shape[0]

  #The Matrix of our laplacian
  if (is_lc):
    B = np.zeros((n,n),PolyElement)
  else:
    B = np.zeros((n,n),float)

  #Sum of adjacent edge weights
  b_x_ys = []
  for i in range(n):
    b_x_y = 0

    for j in range(n):     
        b_x_y = b_x_y + A[i,j]
    b_x_ys.append(b_x_y)
    
    for j in range(n):
      if (i == j):
        B[i,j] = 1
      else:
        #b_x_ys includes A[x,x] right? 
        if (b_x_ys[i] != 0):
            if (isinstance(b_x_ys[i],PolyElement)):
                B[i,j] = -A[i,j] * _series_inversion1((b_x_ys[i]),eps,n_s)
            else:
                B[i,j] = -A[i,j] / (b_x_ys[i])
              
  return B  

def create_dot(A_s : np.array,name_str : str,filename : str):
  #creates a dot file of our graph

  #A_s: adjacency matrix of our graph as a string
  #name_str: name of our file, without file extension
  #create the graph and put the levicvita weights as strings on the edges
  G_l = nx.Graph()

  dim = A_s.shape[0]
  for i in range(0,dim):
    for j in range(0,dim):
      if (A_s[i,j] != '0'):
        G_l.add_edge(i,j,label=A_s[i,j])

  nx.nx_pydot.write_dot(G_l, filename + name_str + ".dot")

def create_graph_example(dim : int, n_eps: int,  max_n: int,n_s : int):
  """
  Creates a graph example with given parameters
  
  Parameters
  ----------
  dim : int
      dimension of the matrix/graph
  n_eps : int
      maximum power of epsilon
  max_n : int
      maximum number for random coefficients

  Creates 6 csv files with and 1 dot file with the lc graph in the directory ./data/
  for example for dim 6, n_eps 2, it creates files of the form 6x6_2_TYPEOFFILE.csv
  in the directory ./data/6x6/eps_2/
  Notes
  -----
  returns an error when the term in the fraction is only a single epsilon
  """

  A,R = create_random_lc_and_real(dim,n_eps,max_n)

  #Calculate the Laplacian
  L = laplacian(A,True,n_s)
  L_r = laplacian(R,False,n_s)
  v = np.linalg.eig(L_r)

  A_s = A.astype(str)
  L_s = L.astype(str)

  dim_str = str(dim)
  dim_x = dim_str + "x" + dim_str
  name_str = dim_x + "_" + str(n_eps)
  filename = "./data/" + dim_x + "/" + "eps_" + str(n_eps) + "/"
  os.makedirs(os.path.dirname(filename), exist_ok=True)
  create_dot(A_s, name_str,filename)

  np.savetxt(filename + name_str + "_adjacency.csv", A_s, delimiter=",",fmt="%s")
  np.savetxt(filename + name_str + "_laplacian.csv", L_s, delimiter=",",fmt="%s")
  np.savetxt(filename + name_str + "_real_adj.csv", R, delimiter=",")
  np.savetxt(filename + name_str + "_real_lap.csv", L_r, delimiter=",")
  np.savetxt(filename + name_str + "_real_eig_vals.csv", v[0], delimiter=",")
  np.savetxt(filename + name_str + "_real_eig_vecs.csv", v[1], delimiter=",")

def create_multiple(max_dim : int, max_eps : int,max_n : int,n_s : int):
  #this creates multiple graphs up to the given dimension max_dim
  #and up to given eps max_eps
  #starting at 3 and 1 resp.

  
  for i in range(2,max_dim):
     for j in range(1,max_eps):
        create_graph_example(i,j,max_n,n_s)