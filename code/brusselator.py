import numpy as np
import matplotlib.pyplot as plt

def brusselator(xl, xr, yb, yt, tb, te, M, N, tsteps):
    """ Program 8.9 Backward differen method for Brusselator
        Input:  xl, xr, yb, yt space interval
                tb, te time interval
                M, N number of space steps in x and y directions
                tsteps number of time steps
        Output: x, y, w mesh and solution """
    Dp, Dq, C, K = 1, 8, 4.5, 7
    delt = (te - tb) / tsteps
    m, n = M + 1, N + 1
    mn = m * n
    mn2 = 2 * mn
    h = (xr - xl)/M
    k = (yt - yb)/N
    x = np.linspace(xl, xr, m)
    y = np.linspace(yb, yt, n)
    fig = plt.figure() 
    ax = fig.add_subplot(111)
    p = np.zeros((m, n))
    q = np.zeros((m, n))
    for i in range(m):          # Define initial conditions
        for j in range(n):
            p[i, j] = C + 0.1
            q[i, j] = K / C + 0.2
    p += 0.*np.random.normal(size=np.shape(p))
    for tstep in range(tsteps):
        v = np.concatenate((p.flatten(), q.flatten()))
        pold = p.copy()
        qold = q.copy()
        for it in range(3):  # Newton iterations
            DF1 = np.zeros((mn2, mn2))
            DF3 = np.zeros((mn2, mn2))
            b = np.zeros(mn2)
            for i in range(1, m-1):
                for j in range(1, n-1):
                    idx = i + j*m
                    DF1[idx, i-1 + j*m] = -Dp/h**2
                    DF1[idx, idx] = Dp*(2/h**2 + 2/k**2) + K+1 + 1/delt
                    DF1[idx, i+1 + j*m] = -Dp/h**2
                    DF1[idx, i + (j-1)*m] = -Dp/k**2
                    DF1[idx, i + (j+1)*m] = -Dp/k**2
                    b[idx] = -pold[i, j]/delt - C
                    DF1[mn + idx, mn + i - 1 + j*m] = -Dq/h**2
                    DF1[mn + idx, mn + idx] = Dq*(2/h**2 + 2/k**2) + 1/delt
                    DF1[mn + idx, mn + i + 1 + j*m] = -Dq/h**2
                    DF1[mn + idx, mn + i + (j - 1)*m] = -Dq/k**2
                    DF1[mn + idx, mn + i + (j+1)*m] = -Dq/k**2
                    DF1[mn + idx, idx] = -K
                    b[mn + idx] = -qold[i, j]/delt
                    DF3[idx, idx] = -2*p[i, j]*q[i, j]
                    DF3[idx, mn + idx] = -p[i, j]**2
                    DF3[mn + idx, idx] = 2*p[i, j]*q[i, j]
                    DF3[mn + idx, mn + idx] = p[i, j]**2
            # Boundary conditions for bottom and top rows
            for i in range(m):
                j = 0 # Bottom
                idx = i + j*m
                DF1[idx, i + j*m] = 3.
                DF1[idx, i + (j+1)*m] = -4.
                DF1[idx, i + (j+2)*m] = 1. 
                DF1[mn+idx, mn+i + j*m] = 3.
                DF1[mn+idx, mn+i + (j+1)*m] = -4.
                DF1[mn+idx, mn+i + (j+2)*m] = 1. 
                j = n-1  # Top
                idx = i + j*m
                DF1[idx, idx] = 3.
                DF1[idx, i + (j-1)*m] = -4.
                DF1[idx, i + (j-2)*m] = 1. 
                DF1[mn+idx, mn+idx] = 3.
                DF1[mn+idx, mn+i + (j-1)*m] = -4.
                DF1[mn+idx, mn+i + (j-2)*m] = 1. 
            # Boundary conditions for left and right columns
            for j in range(1,n-1):
                i = 0 # Left side
                idx = i + j*m
                DF1[idx, idx] = 3.
                DF1[idx, i+1 + j*m] = -4.
                DF1[idx, i+2 + j*m] = 1. 
                DF1[mn+idx, mn+idx] = 3.
                DF1[mn+idx, mn+i+1 + j*m] = -4.
                DF1[mn+idx, mn+i+2 + j*m] = 1. 
                i = m-1  # Right side
                idx = i + j*m
                DF1[idx, idx] = 3.
                DF1[idx, i-1 + j*m] = -4.
                DF1[idx, i-2 + j*m] = 1. 
                DF1[mn+idx, mn+idx] = 3.
                DF1[mn+idx, mn+i-1 + j*m] = -4.
                DF1[mn+idx, mn+i-2 + j*m] = 1. 
            DF = DF1 + DF3   # Combine Jacobian matrices
            F = (DF1 + DF3/3).dot(v) + b
            delta_v = linalg.solve(DF, F)
            v -= delta_v  # Newton step
            p = v[:mn].reshape(m, n)
            q = v[mn:].reshape(m, n)
        X, Y = np.meshgrid(x, y, indexing='ij')  # Plot current solution
        fig = plt.figure(1) 
        fig.clf()
        ax = fig.add_subplot(111)
        ax.contour(X, Y, p, cmap='Blues',levels = 25)
        plt.title('Solution of Brusselator')
        plt.xlim(xl, xr)
        plt.ylim(yb, yt)
        plt.draw()
        plt.pause(.1)
    plt.show()
    return p,q,x,y

# Example usage:
p=brusselator(0,40,0,40,0,10,80,80,100)
