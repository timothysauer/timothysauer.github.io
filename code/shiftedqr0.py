import numpy as np
from scipy.linalg import qr

def shiftedqr0(A):
    """ Program 12.6 Shifted QR Algorithm, preliminary version
        Input:  A square matrix
        Output: lam array of eigenvalues """
    tol = 1e-14
    m = A.shape[0]
    lam = np.zeros(m)  # Initialize eigenvalue array
    n = m
    while n > 1:
        while np.max(np.abs(A[n - 1, :n - 1])) > tol:  # Check 
            mu = A[n - 1, n - 1]  # Define shift mu
            q, r = qr(A - mu*np.eye(n))  # QR factorization
            A = r@q + mu*np.eye(n)  # Update a using QR factors
        lam[n - 1] = A[n - 1, n - 1]  # Declare eigenvalue
        n -= 1  # Decrement n
        A = A[:n, :n]  # Deflate the matrix
    lam[0] = A[0, 0]  # Only a 1x1 matrix remains
    return lam

# Example usage
A = np.array([[16., -6., -15.],   # A muat be symmetric!
              [4., 0., -4.], 
              [10., -4., -9.]])
lam = shiftedqr0(A)
print("Eigenvalues:", lam)
