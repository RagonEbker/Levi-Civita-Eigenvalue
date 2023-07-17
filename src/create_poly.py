from sympy import Poly, symbols
import numpy as np

#This creates the companion Matrix a polynomial with degree 21 with the
#largest root (and therefore eigenvalue of the companion matrix):
# 9*eps**9 + 8*eps**8 + 7*eps**7 + 6*eps**6 + 5*eps**5 + 4*eps**4 + 3*eps**3 + 2*eps**2 + eps + 100
eps = symbols('eps')
lamda = symbols('lamda')
p = Poly((lamda - (100)),lamda)
for i in range(10):
  p -= i*eps**i
degree = 20
for i in range(degree):
  p *= (lamda - (i*2)+(i/degree)*eps**i)

coeffs = np.array(list(reversed(p.coeffs())))
A = np.polynomial.polynomial.polycompanion(coeffs)