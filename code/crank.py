import numpy as np
import matplotlib.pyplot as plt

def crank(xl, xr, yb, yt, M, N):
    """ Program 8.4 Crank-Nicolson Method with
            Dirichlet boundary conditions
        Input:  [xl,xr] space interval
                [yb,yt] time interval
                M number of space steps
                N number of time steps
        Output: w solution """
    def f(x): return np.sin(2*np.pi*x)**2  # Initial condition function
    def l(t): return 0*t  # Left boundary condition
    def r(t): return 0*t  # Right boundary condition
    D = 1  # diffusion coefficient
    h = (xr - xl)/M  # space step size
    k = (yt - yb)/N  # time step size
    sigma = D*k/(h**2)
    m = M - 1  # Number of interior nodes
    n = N  # Number of time steps
    # Define tridiagonal matrices a and b
    a = np.diag(2+2*sigma*np.ones(m))+np.diag(-sigma*np.ones(m-1),1)
    a +=  np.diag(-sigma * np.ones(m - 1), -1)
    b = np.diag(2 - 2*sigma*np.ones(m))+ np.diag(sigma*np.ones(m-1),1) 
    b +=  np.diag(sigma * np.ones(m - 1), -1)
    lside = l(yb + np.arange(n + 1) * k)  # Left boundary values
    rside = r(yb + np.arange(n + 1) * k)  # Right boundary values
    w = np.zeros((m + 2, n + 1))  # Initial conditions
    w[:, 0] = f(xl + np.arange(1, m + 3) * h)
    for j in range(n):   # Time-stepping
        sides = np.zeros(m)
        sides[0] = lside[j] + lside[j + 1]
        sides[1:m - 1] = 0  # Middle part remains zero
        sides[m - 1] = rside[j] + rside[j + 1]
        w[1:m+1, j+1] = np.linalg.solve(a,b@w[1:m+1,j]+sigma*sides)
    w = np.vstack((lside, w, rside))
    x = xl + np.arange(0, M + 3) * h    # grids for plotting
    t = yb + np.arange(0, N + 1) * k
    X, T = np.meshgrid(x, t)
    fig = plt.figure()   # Create a mesh plot
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, T, w.T, cmap='viridis', vmin=-1.,vmax=1.)
    plt.xlabel('x')
    plt.ylabel('t')
    plt.title('Crank-Nicolson Method for heat equation')
    plt.xlim(xl, xr)
    plt.ylim(yb, yt)
    ax.set_zlim(-1.,1.)
    ax.view_init(elev=20, azim=-30)  # Set the view angle
    plt.show()
    return w

# Example usage
w = crank(0, 1, 0, 1, 10, 10)
