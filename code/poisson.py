import numpy as np
import matplotlib.pyplot as plt

def poisson(xl, xr, yb, yt, M, N):
    f = lambda x, y: 0  # Input function data
    g1 = lambda x: np.log(x**2 + 1)  # Boundary values at y = yb
    g2 = lambda x: np.log(x**2 + 4)  # Boundary values at y = yt
    g3 = lambda y: 2 * np.log(y)      # Boundary values at x = xl
    g4 = lambda y: np.log(y**2 + 1)   # Boundary values at x = xr
    m = M + 1
    n = N + 1
    mn = m * n
    h = (xr - xl) / M  # Step size in x
    k = (yt - yb) / N  # Step size in y
    h2 = h ** 2
    k2 = k ** 2
    x = xl + np.arange(0, M + 1) * h  # Mesh values in x
    y = yb + np.arange(0, N + 1) * k  # Mesh values in y
    # Initialize matrix A and vector b
    A = np.zeros((mn, mn))
    b = np.zeros(mn)
    # Fill in the coefficients for interior points
    for i in range(1, m - 1):
        for j in range(1, n - 1):
            A[i + j * m, i - 1 + j * m] = 1 / h2
            A[i + j * m, i + 1 + j * m] = 1 / h2
            A[i + j * m, i + j * m] = -2 / h2 - 2 / k2
            A[i + j * m, i + (j - 1) * m] = 1 / k2
            A[i + j * m, i + (j + 1) * m] = 1 / k2
            b[i + j * m] = f(x[i], y[j])
    # Fill in the coefficients for bottom and top boundary points
    for i in range(m):
        j = 0  # Bottom boundary
        A[i + j * m, i + j * m] = 1
        b[i + j * m] = g1(x[i])
        j = N  # Top boundary
        A[i + j * m, i + j * m] = 1
        b[i + j * m] = g2(x[i])

    # Fill in the coefficients for left and right boundary points
    for j in range(1, n - 1):
        i = 0  # Left boundary
        A[i + j * m, i + j * m] = 1
        b[i + j * m] = g3(y[j])
        
        i = M  # Right boundary
        A[i + j * m, i + j * m] = 1
        b[i + j * m] = g4(y[j])

    # Solve for the solution vector v
    v = np.linalg.solve(A, b)  # Solve Ax = b
    print(v)
    w = v.reshape(m, n)  # Translate from v to w

    # Plotting the solution
    X, Y = np.meshgrid(x, y)
    fig = plt.figure()   # Create a mesh plot
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, w, cmap='viridis', vmin=-1.,vmax=1.)
    plt.xlabel('X')
    plt.ylabel('Y')
    ax.set_zlim(0.,2.)
    plt.title('2D Poisson Equation Solution')
    ax.view_init(elev=20, azim=-135)  # Set the view angle
    plt.show()

    return w

# Example usage
w = poisson(0, 1, 1, 2, 4, 4)
