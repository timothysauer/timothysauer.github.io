import numpy as np

def rqi(A, x, k):
    """ Program 12.3 Rayleigh Quotient Iteration
        Input:  A square matrix
                x initial vector
                k number of steps
        Output: lam eigenvalue
                u eigenvector """
    lam_prev = np.inf
    for j in range(k):
        u = x/np.linalg.norm(x)  # Normalize the vector
        lam = np.dot(u, A@u)      # Rayleigh quotient
        print(lam)
        if np.abs(lam_prev - lam) < 1e-14: break
        lam_prev = lam
        x = np.linalg.solve(A - lam*np.eye(A.shape[0]), u)  
    u = x/np.linalg.norm(x)  # Final normalization
    lam = np.dot(u, A@u)      # Final Rayleigh quotient
    return lam, u

# Example usage
A = np.array([[1., 3.],
              [2., 2.]])
x_initial = np.array([1, 1], dtype=float)  # Initial guess
k_steps = 10  # Number of iterations
lam, u = rqi(A, np.array([1.,3.]), 10)
print("Eigenvalue:", lam)
print("Eigenvector:", u)
