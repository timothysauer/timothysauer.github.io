import numpy as np

def rqi(A, x, k):
    for j in range(k):
        u = x / np.linalg.norm(x)  # Normalize the vector
        lam = np.dot(u, A @ u)      # Rayleigh quotient
        x = np.linalg.solve(A - lam * np.eye(A.shape[0]), u)  # Inverse power iteration

    u = x / np.linalg.norm(x)  # Final normalization
    lam = np.dot(u, A @ u)      # Final Rayleigh quotient
    return lam, u

# Example usage
A = np.array([[4, 1],
              [1, 3]], dtype=float)

x_initial = np.array([1, 1], dtype=float)  # Initial guess
k_steps = 100  # Number of iterations

lambda_value, eigenvector = rqi(A, x_initial, k_steps)
print("Eigenvalue:", lambda_value)
print("Eigenvector:", eigenvector)
