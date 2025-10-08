import numpy as np
import matplotlib.pyplot as plt

def poissonfem(xl, xr, yb, yt, M, N):
    """ Program 8.5 Finite Element Solver for 2d Poisson
        Input:  xl, xr, yb, yt rectangle domain
                M, N space steps
        Output: w solution """
    def f(x): return 0  # Input function of vector x
    def r(x): return 0
    def g1(x): return np.log(x**2 + 1)  # Boundary values at y = yb
    def g2(x): return np.log(x**2 + 4)  # Boundary values at y = yt
    def g3(y): return  2 * np.log(y)      # Boundary values at x = xl
    def g4(y): return np.log(y**2 + 1)   # Boundary values at x = xr
    m, n = M + 1, N + 1
    mn = m * n
    h = (xr - xl) / M
    k = (yt - yb) / N
    h2, k2 = h**2, k**2
    hk = h*k
    x = np.linspace(xl, xr, m)
    y = np.linspace(yb, yt, n)
    A = np.zeros((mn, mn))
    b = np.zeros(mn)
    for i in range(1, m - 1):  # interior points
        for j in range(1, n - 1):
            xi = x[i]
            yj = y[j]
            B1 = [xi - 2*h/3, yj - k/3]
            B2 = [xi - h/3, yj - 2*k/3]
            B3 = [xi + h/3, yj - k/3]
            B4 = [xi + 2*h/3, yj + k/3]
            B5 = [xi + h/3, yj + 2*k/3]
            B6 = [xi - h/3, yj + k/3]
            row = i + j*m
            A[row,row]=2*(h2+k2)/hk-hk*(r(B1)+r(B2)+r(B3)+r(B4)+r(B5)+r(B6))/18
            # Off-diagonal entries
            A[row,(i-1)+j*m] = -k/h - hk*(r(B6)+r(B1))/18
            A[row, (i-1) + (j-1)*m] = -hk*(r(B1)+r(B2))/18
            A[row, i + (j-1)*m] = -hk - hk*(r(B2)+r(B3))/18
            A[row, (i+1)+j*m] = -k/h - hk*(r(B3)+r(B4))/18
            A[row, (i+1)+(j+1)*m] = - hk*(r(B4)*r(B5))/18
            A[row, i+(j+1)*m] = -h/k - hk*(r(B5)+r(B6))/1
            b[row] = -hk*(f(B1)+f(B2)+f(B3)+f(B4)+f(B5)+f(B6))
    for i in range(m):
        j = 0
        row = i+j*m
        A[row,row] = 1
        b[row] = g1(x[i])
        j = n - 1
        row = i+j*m
        A[row,row] = 1
        b[row] = g2(x[i])
    for j in range(1,n-1):
        i = 0
        row = i+j*m
        A[row,row] = 1
        b[row] = g3(y[j])
        i = m-1
        row = i+j*m
        A[row,row] = 1
        b[row] = g4(y[j])
        # Solve for the solution vector v
    v = np.linalg.solve(A, b)  # Solve Ax = b
    w = v.reshape(m, n)  # Translate from v to w
    X, Y = np.meshgrid(x, y)
    fig = plt.figure()   # Create a mesh plot
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, w, cmap='viridis', vmin=-1.,vmax=1.)
    plt.xlabel('X')
    plt.ylabel('Y')
    ax.set_zlim(0.,2.)
    plt.title('2D Poisson Equation Solution by FEM')
    ax.view_init(elev=20, azim=-135)  # Set the view angle
    plt.show()
    return w

# Example usage
w = poissonfem(0, 1, 1, 2, 4, 4)