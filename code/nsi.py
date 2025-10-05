import numpy as np
from scipy.linalg import qr

def nsi(A, k):
    m, n = A.shape
    Q = np.eye(m)  # Initialize Q as an identity matrix
    for j in range(k):
        # Perform QR factorization
        Q, R = qr(A @ Q)
    # Compute eigenvalues using the Rayleigh quotient
    lam = np.diag(Q.T @ A @ Q)
    return lam, Q

# Example usage
A = np.array([[6, 3],
              [3, 4]], dtype=float)
k = 20  # Number of iterations
eigenvalues, eigenvectors = nsi(A, k)
print("Eigenvalues:", eigenvalues)
print("Eigenvector matrix Q:\n", eigenvectors)
