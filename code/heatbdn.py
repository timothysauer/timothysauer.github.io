import numpy as np
import matplotlib.pyplot as plt

def heatbdn(xl, xr, yb, yt, M, N):
    f = lambda x: np.sin(2*np.pi*x)**2  # Initial condition function
    D = 1  # Diffusion coefficient
    h = (xr - xl)/M  # Space step size
    k = (yt - yb)/N  # Time step size
    m = M + 1  # Number of spatial steps
    n = N  # Number of time steps
    sigma = D*k/(h**2)  # Stability parameter
    # Define matrix A
    a = np.diag(1 + 2 * sigma * np.ones(m)) + np.diag(-sigma * np.ones(m - 1), 1) + np.diag(-sigma * np.ones(m - 1), -1)
    a[0,:] = np.concatenate(([-3.,4.,-1.],np.zeros(m-3)))
    a[m-1,:] = np.concatenate((np.zeros(m-3),[-1.,4.,-3]))    # Set boundary conditions
    # Initialize solution matrix
    w = np.zeros((m, n + 1))
    w[:,0] = f(xl + (1+np.arange(m))*h)  # Initial conditions
    # Time-stepping
    for j in range(n):
        w[:,j+1] = np.linalg.solve(a,np.concatenate(([0],w[1:m-1,j],[0])))
    # Attach boundary conditions
    #w = np.vstack((lside, w, rside))
    # Prepare for plotting
    x = np.linspace(xl, xr, m)
    t = np.linspace(yb, yt, n+1)
    # 3-D Plot of the solution
    X, T = np.meshgrid(x, t)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, T, w.T, cmap='viridis', vmin=-1.,vmax=1.)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('w')
    ax.set_title('Heat Equation Solution (Neumann B.C.)')
    plt.xlim(xl, xr)
    plt.ylim(yb, yt)
    ax.set_zlim(-1.,1.)
    ax.view_init(elev=30, azim=-30)  # Set the view angle
    plt.show()

    return w

# Example usage
w = heatbdn(0, 1, 0, 1, 20, 20)

