from sympy.polys.ring_series import *
from sympy import *
import numpy as np

R, eps = ring('eps', RR)

def str_to_lc(lc_s : str):
        lc = parse_expr(lc_s,local_dict={'eps':eps})
        return lc

def read_lc_matrix(filename : str):
    csv = np.array(np.genfromtxt(filename,delimiter=',',dtype=str).tolist())
    lc_array = np.array([[str_to_lc(el) for el in row] for row in csv])
    return lc_array
