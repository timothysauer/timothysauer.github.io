import numpy as np
import matplotlib.pyplot as plt

def euler2(ydot,inter, y0, n):
    """Program 6.2 Euler's Method, vector version
    Solve the IVP y' = ydot(t,y)
    Input:  ydot right-hand side function of ODE
            inter time interval of solution
            y0 initial vector
            n number of time steps
    Output: t time points
            y solution"""
    t = np.zeros(n + 1)  # Time steps
    y = np.zeros((n + 1, len(y0)))  # Solution 
    t[0] = inter[0]
    y[0, :] = y0
    h = (inter[1] - inter[0]) / n  # Step size
    for i in range(n):
        t[i + 1] = t[i] + h
        y[i + 1, :] = eulerstep(t[i], y[i, :], h)
    plt.plot(t, y[:, 0], label='y0')
    plt.plot(t, y[:, 1], label='y1')
    plt.xlabel('Time')
    plt.ylabel('Solution')
    plt.title('Euler Method')
    plt.legend()
    plt.show()
    return t, y

def eulerstep(t, y, h):
    return y + h * ydot(t, y)

# Example usage
#def ydot(t, y):
#    z = np.zeros(2)
#    z[0] = y[1]**2 - 2 * y[0]
#    z[1] = y[0] - y[1] - t * y[1]**2
#    return z
#y = euler2(ydot,[0, 1], [0, 1], 10)
