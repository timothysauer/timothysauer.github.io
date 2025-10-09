import numpy as np

def powerit(A, x, k):
    """ Program 12.1 Power Iteration
        Computes domainant eignvalue/eigenvector
        Input:  A matrix
                x initial vector
                n number of steps
        Output: lam dominant eigenvalue
                u domeinant eigenvector """
    for j in range(k):
        u = x/np.linalg.norm(x)  # Normalize the vector
        x = A@u  # Power step
        lam = u.T@x  # Rayleigh Quotient
    u = x/np.linalg.norm(x)  # Final normalization 
    return lam, u

# Example usage
A = np.array([[1., 3.],
              [2., 2.]])
lam, u = powerit(A,np.array([1.,1.]),10)

