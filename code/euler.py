import numpy as np
import matplotlib.pyplot as plt

def euler(inter, y0, n):
    t = np.zeros(n + 1)  # Time steps
    y = np.zeros(n + 1)  # Solution 
    t[0] = inter[0]
    y[0] = y0
    h = (inter[1] - inter[0]) / n  # Step size

    for i in range(n):
        t[i + 1] = t[i] + h
        y[i + 1] = eulerstep(t[i], y[i], h)

    # Plotting the results
    plt.plot(t, y)
    plt.xlabel('Time')
    plt.ylabel('Solution')
    plt.title('Euler Method')  
    plt.show()

    return t, y

def eulerstep(t, y, h):
    """ One step of the Euler method. """
    return y + h * ydot(t, y)

def ydot(t, y):
    """ Right-hand side of differential equation """ 
    return t*y + t**3

# Example usage
y = euler([0, 1], 0.1, 10)
