import numpy as np
from scipy.sparse import diags

def sparsesetup(n):
    e = np.ones(n)
    n2 = n // 2
    diagonals = [-e, 3*e, -e] # Create the sparse matrix a using diagonals
    offsets = [-1, 0, 1]
    a = diags(diagonals, offsets, shape=(n, n)).tocsc() 
    # Create the left half diagonal and add 
    c = diags([e / 2], [0], shape=(n, n)).tocsc()
    a += c[:,::-1]  # Flip the left half diagonal to the right
    a[n2, n2-1] = -1       # Fix up two entries
    a[n2-1, n2] = -1
    b = np.zeros(n)  # Create the right-hand side vector 'b'
    b[0] = 2.5
    b[n-1] = 2.5
    b[1:n-1] = 1.5
    b[n2-1:n2+1] = 1.

    return a, b

