import numpy as np

def hessen(a):
    m, n = a.shape
    v = np.zeros((m, m))
    for k in range(m - 2):
        x = a[k+1:m+1,k]  # Vector for Householder transformation
        x_sign = -np.sign(x[0]+np.finfo(float).eps)
        v[0:m-k-1,k] = x_sign*np.linalg.norm(x)*np.eye(m-k-1,1).squeeze() - x
        v[0:m-k-1,k] = v[0:m-k-1,k]/np.linalg.norm(v[0:m-k-1, k])  # Normalize v
        # Apply reflection to a, on both sides
        a[k+1:m,k:m] -= 2*np.outer(v[0:m-k-1,k],v[0:m-k-1,k]) @ a[k+1:m,k:m]
        a[:,k+1:m] -= 2*a[:,k+1:m] @ np.outer(v[0:m-k-1,k],v[0:m-k-1, k])
    return a, v

# Example usage
A = np.array([[4, 1, 2,3], 
              [1, 2, 3, 4], 
              [2, 3, 6, 3],
              [3, 1, -3, 4]], dtype=float)
H, V = hessen(A)
print("Hessenberg form:\n", H)
print("Reflectors:\n", V)
