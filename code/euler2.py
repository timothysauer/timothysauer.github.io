import numpy as np
import matplotlib.pyplot as plt

def euler2(inter, y0, n):
    t = np.zeros(n + 1)  # Time steps
    y = np.zeros((n + 1, len(y0)))  # Solution 
    t[0] = inter[0]
    y[0, :] = y0
    h = (inter[1] - inter[0]) / n  # Step size

    for i in range(n):
        t[i + 1] = t[i] + h
        y[i + 1, :] = eulerstep(t[i], y[i, :], h)

    # Plotting the results
    plt.plot(t, y[:, 0], label='y0')
    plt.plot(t, y[:, 1], label='y1')
    plt.xlabel('Time')
    plt.ylabel('Solution')
    plt.title('Euler Method')
    plt.legend()
    plt.show()

    return t, y

def eulerstep(t, y, h):
    """ One step of the Euler method. """
    return y + h * ydot(t, y)

def ydot(t, y):
    """ Function that defines the derivatives. """
    z = np.zeros(2)
    z[0] = y[1]**2 - 2 * y[0]
    z[1] = y[0] - y[1] - t * y[1]**2
    return z

# Example usage
y = euler2([0, 1], [0, 1], 10)
