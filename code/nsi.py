import numpy as np
from scipy.linalg import qr

def nsi(A, k):
    """ Program 12.4 Normlized Simultaneous Iteration
        Input:  A square matrix
                k number of steps
        Output: lam array of eigenvalues
                u matrix of eigenvectors as columns """
    m, n = A.shape
    Q = np.eye(m)  # Initialize Q as an identity matrix
    for j in range(k):
        # Perform QR factorization
        Q, R = qr(A@Q)
    # Compute eigenvalues using the Rayleigh quotient
    lam = np.diag(Q.T@A@Q)
    return lam, Q

# Example usage
A = np.array([[1., 3.],
              [2., 2.]], dtype=float)
eigenvalues, eigenvectors = nsi(A, 20)
print("Eigenvalues:", eigenvalues)
print("Eigenvector matrix Q:\n", eigenvectors)
