import numpy as np
from scipy.linalg import qr

def shifted_qr(a):
    tol = 1e-14
    kounttol = 500
    m = a.shape[0]
    lam = np.zeros(m,dtype='complex')  # Initialize complex eigenvalue array
    n = m
    while n > 1:
        kount = 0
        while np.max(np.abs(a[n-1,:n-1])) > tol and kount < kounttol:
            kount += 1  # Keep track of the number of QR iterations
            mu = a[n-1,n-1]  # Shift is mu
            q, r = qr(a - mu*np.eye(n))  # QR factorization
            a = r @ q + mu*np.eye(n)  # Update a using the QR factors
        if kount < kounttol:  # Have isolated 1x1 block
            lam[n-1] = a[n-1,n-1]  # Declare eigenvalue
            n -= 1
            a = a[:n, :n]  # Deflate by 1
        else:  # Have isolated 2x2 block
            disc = (a[n-2,n-2] - a[n-1,n-1])**2 + 4*a[n-1,n-2]*a[n-2,n-1]
            lam[n-1] = (a[n-2,n-2] + a[n-1,n-1] + np.emath.sqrt(disc))/2
            lam[n-2] = (a[n-2,n-2] + a[n-1,n-1] - np.emath.sqrt(disc))/2
            n -= 2
            a = a[:n, :n]  # Deflate by 2
    if n > 0:
        lam[0] = a[0, 0]  # Only a 1x1 block remains
    return lam

# Example usage
A = np.array([[3,4,6],
              [5,-1,-1],
              [-15,22,16]])
eigenvalues = shifted_qr(A)
print("Eigenvalues:", eigenvalues)
