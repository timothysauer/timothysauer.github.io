import numpy as np

def houseqr(A):
    """
    Program 4.2 Orthogonalization by Householder reflectors.
    Input: 
        A: mxn matrix with linearly independent columns
    Output: 
        Q: orthogonal matrix
        R: upper triangular matrix
    """
    m, n = A.shape
    Q = np.eye(m)  # Initialize Q as the identity matrix
    for i in range(min(n, m - 1)):
        # Extract the vector x
        x = A[i:m, i]
        # Calculate the Householder vector
        w = np.zeros_like(x)
        w[0] = -np.sign(x[0]) * np.linalg.norm(x)
        v = w - x
        # Form the Householder matrix
        H = np.eye(m)
        H[i:m, i:m] = np.eye(m - i) - 2 * np.outer(v, v) / np.dot(v, v)
        # Update Q and A
        Q = Q @ H  # Q is updated as the product of Householder transformations
        A = H @ A  # Apply the Householder transformation to A

    R = A  # R is the upper triangular part of the transformed A
    return Q, R
