import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

def bvpfem(inter, bv, n):
    a = inter[0]
    b = inter[1]
    ya = bv[0]
    yb = bv[1]
    h = (b - a) / (n + 1)
    alpha = (8 / 3) * h + 2 / h
    beta = (2 / 3) * h - 1 / h
    # Create the sparse matrix M
    e = np.ones(n)
    diagonals = [beta * e, alpha * e, beta * e]
    M = diags(diagonals, offsets=[-1, 0, 1], shape=(n, n),format='csr')
    # Create the right-hand side vector d
    d = np.zeros(n)
    d[0] = -ya * beta
    d[-1] = -yb * beta
    # Solve the linear system M * c = d
    c = spsolve(M, d)
    # Plotting w with boundary values
    plt.plot([a] + list(a+(1 + np.arange(n))* +[b],[ya] + list(c) + [yb])
    plt.xlabel('x')
    plt.ylabel('w')
    plt.title(' BVP FEM Solution')
    plt.grid(True)
    plt.show()
    return c

# Example usage
c = bvpfem([0, 1], [1, 3], 9)
print("Solution values:", c)
