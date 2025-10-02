import numpy as np

def hessen(a):
    m, n = a.shape
    v = np.zeros((m, m))

    for k in range(m - 2):
        x = a[k+1:m, k]  # Vector for Householder transformation
        v[0:m-k, k] = -np.sign(x[0] + np.finfo(float).eps) * np.linalg.norm(x) * np.eye(m-k, 1).squeeze() - x
        v[0:m-k, k] = v[0:m-k, k] / np.linalg.norm(v[0:m-k, k])  # Normalize v

        # Apply reflection to a
        a[k + 1:m, k:m] -= 2 * np.outer(v[0:m-k, k], v[0:m-k, k].T @ a[k + 1:m, k:m])
        a[0:m, k + 1:m] -= 2 * (a[:, k + 1:m] @ v[0:m-k, k])[:, np.newaxis] * v[0:m-k, k]

    return a, v

# Example usage
A = np.array([[4, 1, 2], 
              [1, 2, 3], 
              [2, 3, 6]], dtype=float)

H, V = hessen(A)
print("Hessenberg form:\n", H)
print("Reflectors:\n", V)
