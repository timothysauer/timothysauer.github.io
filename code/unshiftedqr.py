import numpy as np
from scipy.linalg import qr

def unshifted_qr(A, k):
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
# A = np.array([6,3],[3,8]]) 
# eigvals, eigvecs = unshifted_qr(A, 10)