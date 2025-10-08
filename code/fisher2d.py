import numpy as np
import matplotlib.pyplot as plt

def fisher2d(xl, xr, yb, yt, tb, te, M, N, tsteps):
    """ Program 8.8 Backward Difference Method for Fisher's equation
            on 2D domain
        Input:  xl, xr, yb, yt space interval
                tb, te time interval
                M, N number of space steps in x and y directions
                tsteps number of time steps
        Output: x, y, w mesh and solution """
    def f(x,y): return 2 + np.cos(np.pi * x) * np.cos(np.pi * y)
    delt = (te - tb) / tsteps
    D = 1.0
    m = M + 1
    n = N + 1
    mn = m * n
    h = (xr - xl) / M
    k = (yt - yb) / N
    x = np.linspace(xl, xr, m)
    y = np.linspace(yb, yt, n)
    fig = plt.figure() 
    w = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            w[i, j] = f(x[i], y[j])
    for tstep in range(tsteps):
        v = w.reshape(mn)
        wold = w.copy()
        for it in range(3):  # Newton iterations
            b = np.zeros(mn)
            DF1 = np.zeros((mn, mn))
            DF2 = np.zeros((mn, mn))
            for i in range(1, m - 1):  # Build Jacobian and RHS
                for j in range(1, n - 1):
                    idx = i + j * m
                    DF1[idx, i - 1 + j*m] = -D/h**2
                    DF1[idx, i + 1 + j*m] = -D/h**2
                    DF1[idx, idx] = 2*D/h**2 + 2*D/k**2 - 1 + 1/delt
                    DF1[idx, i + (j - 1)*m] = -D/k**2
                    DF1[idx, i + (j + 1)*m] = -D/k**2
                    b[idx] = -wold[i, j]/delt
                    DF2[idx,idx] = 2*w[i, j]
            # Boundary conditions for bottom and top rows
            for i in range(m):
                j = 0 # Bottom
                idx = i + j*m
                DF1[idx, i + j*m] = 3
                DF1[idx, i + (j+1)*m] = -4
                DF1[idx, i + (j+2)*m] = 1 
                j = n-1  # Top
                idx = i + j*m
                DF1[idx, idx] = 3
                DF1[idx, i + (j-1)*m] = -4
                DF1[idx, i + (j-2)*m] = 1 
            # Boundary conditions for left and right columns
            for j in range(1,n-1):
                i = 0 # Left side
                idx = i + j*m
                DF1[idx, idx] = 3
                DF1[idx, i+1 + j*m] = -4
                DF1[idx, i+2 + j*m] = 1 
                i = m-1  # Right side
                idx = i + j*m
                DF1[idx, idx] = 3
                DF1[idx, i-1 + j*m] = -4
                DF1[idx, i-2 + j*m] = 1 
            DF = DF1 + DF2 # Combine Jacobian matrices
            F = (DF1 + DF2/2).dot(v) + b
            delta_v = np.linalg.solve(DF, F) # Newton step
            v -= delta_v
            w = v.reshape(m, n)
        X, Y = np.meshgrid(x, y, indexing='ij')
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, Y, w.T, cmap='Blues', vmin=-1.,vmax=1.)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Solution of Fisher 2D Equation')
        plt.xlim(xl, xr)
        plt.ylim(yb, yt)
        ax.set_zlim(0.,3)
        ax.view_init(elev=20, azim=-120)  # Set the view angle
        plt.draw()
        plt.pause(.05)
    plt.show()
    return X, Y, w

# Example usage
w = fisher2d(0,1,0,1,0,1,20,20,25)
