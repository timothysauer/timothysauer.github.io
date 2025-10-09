import numpy as np
from scipy.linalg import qr

def unshifted_qr(A, k):
    """ Program 12.5 Unshifted QR Algorithm
        Input:  A square matrix
                k number of steps
        Output: lam array of eigenvalues
                Qbar matrix of eigenvalues as columns """
    m, n = A.shape
    Q = np.eye(m)
    Qbar = Q.copy()
    R = A.copy()
    for _ in range(k):
        Q, R = qr(R @ Q)
        Qbar = Qbar @ Q
    lam = np.diag(R @ Q)
    return lam, Qbar

# Example usage:
A = np.array([[1.,3.],[2.,2.]]) 
lam, Qbar = unshifted_qr(A, 20)
print("Eigenvalue:", lam)
print("Eigenvector:", Qbar)