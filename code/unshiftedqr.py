import numpy as np

def unshifted_qr(A, k):
    m, n = A.shape
    Q = np.eye(m)
    Qbar = Q.copy()
    R = A.copy()

    for _ in range(k):
        Q, R = np.linalg.qr(R @ Q)
        Qbar = Qbar @ Q

    lam = np.diag(R @ Q)
    return lam, Qbar

# Example usage:
# A = np.array([[...], [...], ...])  # your symmetric matrix
# eigvals, eigvecs = unshifted_qr(A, 100)