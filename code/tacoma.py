import numpy as np
import matplotlib.pyplot as plt

def tacoma(inter, ic, n, p):
    plt.ion()  # Enable interactive plotting
    plt.figure(figsize=(8, 8))
    plt.axis([-8, 8, -8, 8])
    plt.xlabel('X Position')
    plt.ylabel('Y Position')

    h = (inter[1] - inter[0]) / n
    y = np.zeros((n + 1, len(ic)))  # Initialize solution array
    y[0, :] = ic  # Set initial conditions
    t = np.zeros(n + 1)
    t[0] = inter[0]
    length = 6  # Length of the bridge

    # Prepare lines for the road and cables
    road_line, = plt.plot([], [], 'b-', linewidth=5)
    left_cable_line, = plt.plot([], [], 'r-', linewidth=1)
    right_cable_line, = plt.plot([], [], 'r-', linewidth=1)

    for k in range(n):
        for i in range(p):
            t[i + 1] = t[i] + h
            y[i + 1, :] = trapstep(t[i], y[i, :], h)

        # Update the state for the next visualization step
        y[0, :] = y[p, :]
        t[0] = t[p]
        
        c = length * np.cos(y[0, 2])  # Calculate horizontal cable position
        s = length * np.sin(y[0, 2])  # Calculate vertical cable position
        
        # Update positions of road and cables
        road_line.set_data([-c, c], [-s - y[0, 0], s - y[0, 0]])
        left_cable_line.set_data([-c, -c], [-s - y[0, 0], 8])
        right_cable_line.set_data([c, c], [s - y[0, 0], 8])
        
        plt.draw()
        plt.pause(h)

    plt.ioff()  # Disable interactive mode
    plt.show()

def trapstep(t, x, h):
    """ One step of the Trapezoid Method """
    z1 = ydot(t, x)
    g = x + h * z1
    z2 = ydot(t + h, g)
    return x + h * (z1 + z2) / 2

def ydot(t, y):
    """ Compute the derivatives for the IVP """
    length = 6
    a = 0.2
    W = 80
    omega = 2 * np.pi * 38 / 60
    a1 = np.exp(a * (y[0] - length * np.sin(y[2])))
    a2 = np.exp(a * (y[0] + length * np.sin(y[2])))
    
    dydt = np.zeros(4)
    dydt[0] = y[1]
    dydt[1] = -0.01 * y[1] - 0.4 * (a1 + a2 - 2) / a + 0.2 * W * np.sin(omega * t)
    dydt[2] = y[3]
    dydt[3] = -0.01 * y[3] + 1.2 * np.cos(y[2]) * (a1 - a2) / (length * a)
    
    return dydt

# Example usage
tacoma([0, 100], [1, 0, 0.001, 0], 25000, 5)
