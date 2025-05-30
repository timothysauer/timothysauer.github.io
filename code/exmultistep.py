import numpy as np
import matplotlib.pyplot as plt

def exmultistep(inter, ic, n, s):
    h = (inter[1] - inter[0]) / n
    y = np.zeros((n + 1, len(ic)))  # Initialize solution array
    f = np.zeros((n + 1, len(ic)))  # Store derivative evaluations
    t = np.zeros(n + 1)  # Initialize time steps

    # Start-up phase
    y[0, :] = ic
    t[0] = inter[0]
    
    for i in range(s - 1):  # start-up phase, using the one-step method
        t[i + 1] = t[i] + h
        y[i + 1, :] = trapstep(t[i], y[i, :], h)
        f[i, :] = ydot(t[i], y[i, :])

    for i in range(s - 1, n):  # multistep method loop
        t[i + 1] = t[i] + h
        f[i, :] = ydot(t[i], y[i, :])
        y[i + 1, :] = ab2step(t[i], i, y, f, h)

    # Plot the results
    plt.plot(t, y)
    plt.xlabel('Time')
    plt.ylabel('Solution')
    plt.title('Multistep Method Solution')
    plt.legend(['y1', 'y2'])  # Change according to the number of variables
    plt.show()

    return t, y

def trapstep(t, x, h):
    """ One step of the Trapezoid Method """
    z1 = ydot(t, x)
    g = x + h * z1
    z2 = ydot(t + h, g)
    return x + h * (z1 + z2) / 2

def ab2step(t, i, y, f, h):
    """ One step of the Adams-Bashforth 2-step method """
    return y[i, :] + h * (3 * f[i, :] / 2 - f[i - 1, :] / 2)

def unstable2step(t, i, y, f, h):
    """ One step of an unstable 2-step method """
    return -y[i, :] + 2 * y[i - 1, :] + h * (5 * f[i, :] / 2 + f[i - 1, :] / 2)

def weaklystable2step(t, i, y, f, h):
    """ One step of a weakly-stable 2-step method """
    return y[i - 1, :] + h * 2 * f[i, :]

def ydot(t, y):
    """ Define the derivative function for the IVP """
    return t * y + t**3

# Example usage
t, y = exmultistep([0, 1], [1], 20, 2)
