import numpy as np
from scipy.linalg import qr

def shiftedqr0(a):
    tol = 1e-14
    m = a.shape[0]
    lam = np.zeros(m)  # Initialize eigenvalue array
    n = m
    while n > 1:
        while np.max(np.abs(a[n - 1, :n - 1])) > tol:  # Check 
            mu = a[n - 1, n - 1]  # Define shift mu
            q, r = qr(a - mu * np.eye(n))  # QR factorization
            a = r @ q + mu * np.eye(n)  # Update a using QR factors
        lam[n - 1] = a[n - 1, n - 1]  # Declare eigenvalue
        n -= 1  # Decrement n
        a = a[:n, :n]  # Deflate the matrix
    lam[0] = a[0, 0]  # Only a 1x1 matrix remains

    return lam

# Example usage
A = np.array([[4, 1, 2],   # A muat be symmetric!
              [1, 3, 0], 
              [2, 0, 2]], dtype=float)

eigenvalues = shiftedqr0(A)
print("Eigenvalues:", eigenvalues)
