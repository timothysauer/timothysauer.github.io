import numpy as np

def powerit(A, x, k):
    for j in range(k):
        u = x / np.linalg.norm(x)  # Normalize the vector
        x = A @ u  # Power step
        lam = u.T @ x  # Rayleigh Quotient
    u = x / np.linalg.norm(x)  # Final normalization 
    return lam, u

# Example usage
A = np.array([[1, 3],
              [2, 2]], dtype=float)
x_initial = np.array([1, 1], dtype=float)  # Initial guess
k_steps = 10  # Number of iterations
lambda_value, eigenvector = powerit(A, x_initial, k_steps)
print("Dominant Eigenvalue:", lambda_value)
print("Dominant Eigenvector:", eigenvector)
