import numpy as np
import matplotlib.pyplot as plt

def heatfd(xl, xr, yb, yt, M, N):
    f = lambda x: np.sin(2*np.pi*x)**2  # Initial condition function
    l = lambda t: 0*t  # Left boundary condition
    r = lambda t: 0*t  # Right boundary condition
    D = 1  # Diffusion coefficient

    h = (xr - xl)/M  # Space step size
    k = (yt - yb)/N  # Time step size
    m = M - 1  # Number of spatial steps
    n = N  # Number of time steps
    sigma = D*k/(h**2)  # Stability parameter
    # Define matrix A
    a = np.diag(1-2*sigma*np.ones(m))+np.diag(sigma*np.ones(m-1),1)+np.diag(sigma*np.ones(m-1),-1)
    # Set boundary conditions
    lside = l(yb + np.arange(n+1)*k)  # Left boundary values
    rside = r(yb + np.arange(n+1)*k)  # Right boundary values
    # Initialize solution matrix
    w = np.zeros((m, n + 1))
    w[:,0] = f(xl + (1+np.arange(m))*h)  # Initial conditions
    # Time-stepping
    for j in range(n):
        w[:,j+1] = a@w[:,j] + sigma*np.concatenate(([lside[j]],np.zeros(m-2),[rside[j]]))
    # Attach boundary conditions
    w = np.vstack((lside, w, rside))
    # Prepare for plotting
    x = np.linspace(xl, xr, m + 2)
    t = np.linspace(yb, yt, n + 1)
    # 3-D Plot of the solution
    X, T = np.meshgrid(x, t)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, T, w.T, cmap='viridis', vmin=-1.,vmax=1.)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('w')
    ax.set_title('Heat Equation Solution (Forward Difference Method)')
    plt.xlim(xl, xr)
    plt.ylim(yb, yt)
    ax.set_zlim(-1.,1.)
    ax.view_init(elev=30, azim=-30)  # Set the view angle
    plt.show()

    return w

# Example usage
w = heatfd(0, 1, 0, 1, 10, 250)
