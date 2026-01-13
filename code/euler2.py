import numpy as np
import matplotlib.pyplot as plt

def euler2(f,inter, y0, n):
    """Program 6.2 Euler's Method, vector version
    Solve the IVP y' = f(t,y)
    Input:  f right-hand side function of ODE
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
        y[i + 1, :] = eulerstep(f, t[i], y[i, :], h)
    return t, y

def eulerstep(f, t, y, h): return y + h * f(t, y)

# Example usage
def f(t, y):
    z = np.zeros(2)
    z[0] = y[1]**2 - 2 * y[0]
    z[1] = y[0] - y[1] - t * y[1]**2
    return z
t, y = euler2(f, [0, 1], [0, 1], 10)
plt.plot(t, y[:, 0], label='y0')
plt.plot(t, y[:, 1], label='y1')
plt.xlabel('Time')
plt.ylabel('Solution')
plt.title('Euler Method')
plt.legend()
plt.show()
