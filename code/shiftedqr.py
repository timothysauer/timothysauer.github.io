import numpy as np
from scipy.linalg import qr

def shiftedqr(A):
    """ Program 12.7 Shifted QR Algorithm
        Input:  A matrix
        Output: lam array of eigenvalues """
    tol = 1e-14
    kounttol = 500
    m = A.shape[0]
    lam = np.zeros(m,dtype='complex')  # Initialize complex eigenvalue array
    n = m
    while n > 1:
        kount = 0
        while np.max(np.abs(A[n-1,:n-1])) > tol and kount < kounttol:
            kount += 1  # Keep track of the number of QR iterations
            mu = A[n-1,n-1]  # Shift is mu
            q, r = qr(A - mu*np.eye(n))  # QR factorization
            A = r @ q + mu*np.eye(n)  # Update a using the QR factors
        if kount < kounttol:  # Have isolated 1x1 block
            lam[n-1] = A[n-1,n-1]  # Declare eigenvalue
            n -= 1
            A = A[:n, :n]  # Deflate by 1
        else:  # Have isolated 2x2 block
            disc = (A[n-2,n-2] - A[n-1,n-1])**2 + 4*A[n-1,n-2]*A[n-2,n-1]
            lam[n-1] = (A[n-2,n-2] + A[n-1,n-1] + np.emath.sqrt(disc))/2
            lam[n-2] = (A[n-2,n-2] + A[n-1,n-1] - np.emath.sqrt(disc))/2
            n -= 2
            A = A[:n, :n]  # Deflate by 2
    if n > 0:
        lam[0] = A[0, 0]  # Only a 1x1 block remains
    return lam

# Example usage
A = np.array([[3.,4.,6.],
              [5.,-1.,-1.],
              [-15.,22.,16.]])
lam = shiftedqr(A)
print("Eigenvalues:", lam)
