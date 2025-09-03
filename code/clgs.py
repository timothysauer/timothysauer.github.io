import numpy as np

def clgs(A):
    """
    Classical Gram-Schmidt orthogonalization.
    
    Input: 
        A: mxn matrix with linearly independent columns
    Output: 
        Q: orthogonal matrix
        R: upper triangular matrix
    """
    m, n = A.shape
    Q = np.zeros((m, n))
    R = np.zeros((n, n))
    
    for j in range(n):
        y = A[:, j].copy()
        for i in range(j):
            R[i, j] = np.dot(Q[:, i], A[:, j])  # Inner product
            y -= R[i, j] * Q[:, i]              # Subtract the projection
        R[j, j] = np.linalg.norm(y)            # Norm of y
        Q[:, j] = y / R[j, j]                  # Normalize y to form the j-th column of Q

    return Q, R
