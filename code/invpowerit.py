import numpy as np

def invpowerit(A, x, s, k):
    """ Program 12.2 Inverse Power Iteration
        Computes eigenvalue nearest to shift s
        Input:  A square matrix
                x innitial vector
                s shift
                k number of steps
        Output: lam eigenvalue closest to s
                u eigenvector """
    As = A - s*np.eye(A.shape[0])  # Shifted matrix A - s*I
    for j in range(k):
        u = x/np.linalg.norm(x)  # Normalize vector
        x = np.linalg.solve(As,u)  # Power step
        lam = np.dot(u, x)  # Rayleigh Quotient
    lam = 1/lam + s  # Adjust for the shift
    u = x/np.linalg.norm(x)  # Normalize the resulting eigenvector
    return lam, u

# Example usage
A = np.array([[1,3],
              [2, 2]], dtype=float)
lam, u = invpowerit(A,np.array([1.,0.]),2.0,20)


