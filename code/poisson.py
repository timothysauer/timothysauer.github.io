import numpy as np
import matplotlib.pyplot as plt

def poisson(xl, xr, yb, yt, M, N):
    """ Program 8.5 Finite Difference Solver for 2d Poisson
        Input:  xl, xr, yb, yt rectangle domain
                M, N space steps
        Output: w solution """
    def f(x,y): return 0  # Input function data
    def g1(x): return np.log(x**2 + 1)  # Boundary values at y = yb
    def g2(x): return np.log(x**2 + 4)  # Boundary values at y = yt
    def g3(y): return  2 * np.log(y)      # Boundary values at x = xl
    def g4(y): return np.log(y**2 + 1)   # Boundary values at x = xr
    m, n = M + 1, N +1
    mn = m*n
    h = (xr - xl)/M  # Step size in x
    k = (yt - yb)/N  # Step size in y
    h2, k2 = h**2, k**2
    x = xl + np.arange(0, M + 1) * h  # Mesh values in x
    y = yb + np.arange(0, N + 1) * k  # Mesh values in y
    A = np.zeros((mn, mn))
    b = np.zeros(mn)
    for i in range(1, m - 1):   # Fill in coefficients for interior points
        for j in range(1, n - 1):
            A[i + j*m, i - 1 + j*m] = 1/h2
            A[i + j*m, i + 1 + j*m] = 1/h2
            A[i + j*m, i + j*m] = -2/h2 - 2/k2
            A[i + j*m, i + (j - 1)*m] = 1/k2
            A[i + j*m, i + (j + 1)*m] = 1/k2
            b[i + j*m] = f(x[i], y[j])
    for i in range(m):
        j = 0  # Bottom boundary
        A[i + j*m, i + j*m] = 1
        b[i + j*m] = g1(x[i])
        j = N  # Top boundary
        A[i + j*m, i + j*m] = 1
        b[i + j*m] = g2(x[i])
    for j in range(1, n - 1):
        i = 0  # Left boundary
        A[i + j*m, i + j*m] = 1
        b[i + j*m] = g3(y[j])
        i = M  # Right boundary
        A[i + j*m, i + j*m] = 1
        b[i + j*m] = g4(y[j])
    v = np.linalg.solve(A, b)  # Solve Ax = b
    w = v.reshape(m, n)  # Translate from v to w
    X, Y = np.meshgrid(x, y)
    return X, Y, w

# Example usage
X, Y, w = poisson(0, 1, 1, 2, 4, 4)
fig = plt.figure()   # Create a mesh plot
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, w, cmap='viridis', vmin=-1.,vmax=1.)
plt.xlabel('X')
plt.ylabel('Y')
ax.set_zlim(0.,2.)
plt.title('2D Poisson Equation Solution')
ax.view_init(elev=20, azim=-135)  # Set the view angle
plt.show()