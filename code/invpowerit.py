import numpy as np

def invpowerit(A, x, s, k):
    As = A - s * np.eye(A.shape[0])  # Shifted matrix A - s*I
    for j in range(k):
        u = x / np.linalg.norm(x)  # Normalize vector
        x = np.linalg.solve(As, u)  # Power step
        lam = np.dot(u, x)  # Rayleigh Quotient
    lam = 1 / lam + s  # Adjust for the shift
    u = x / np.linalg.norm(x)  # Normalize the resulting eigenvector

    return lam, u

# Example usage
A = np.array([[4, 1],
              [2, 3]], dtype=float)
x = np.array([1, 0], dtype=float)  # Initial guess
s = 1.0  # Shift
k = 10   # Number of steps


