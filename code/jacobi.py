import numpy as np
from scipy.sparse import diags,issparse

def jacobi(a, b, k):
    d = a.diagonal()  # Extract diagonal of a
    r = a - diags(d)  # r is the remainder
    if not issparse(a):
        r = np.array(r)  
    x = np.zeros((len(b)))  # Initialize vector x
    for j in range(k):  # Loop for Jacobi iteration
        x = (b - r@x)/d 
    return x

