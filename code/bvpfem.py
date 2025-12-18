import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

def bvpfem(inter, bv, n):
    """ Program 7.2 Finite element solution fo linear BVP
        Input:  inter time interval
                bv boundary values
                n number of steps
        Output: w solution """
    a = inter[0]
    b = inter[1]
    ya = bv[0]
    yb = bv[1]
    h = (b - a)/(n + 1)
    alpha = (8/3)*h + 2/h
    beta = (2/3)*h - 1/h
    # Create the sparse matrix M
    e = np.ones(n)
    diagonals = [beta*e, alpha*e, beta*e]
    M = diags(diagonals, offsets=[-1,0,1], shape=(n, n),format='csr')
    # Create the right-hand side vector d
    d = np.zeros(n)
    d[0] = -ya*beta
    d[-1] = -yb*beta
    w = spsolve(M, d)  # Solve the linear system M * c = d
    t = [a] + list(a+(1+np.arange(n))*h)+[b]
    w = [ya]+list(w)+[yb]
    return t, w

# Example usage
t, y = bvpfem([0, 1], [1, 3], 39)
plt.xlabel('t')
plt.ylabel('y')
plt.title(' BVP Finite Element Method Solution')
plt.grid(True)
plt.show()
