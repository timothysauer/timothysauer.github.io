import numpy as np
import matplotlib.pyplot as plt

def tacoma(inter, ic, n, p):
    """ Program 6.6 Animation of Tacoma Narrows Bridge
        Input:  inter time interval
                ic initial conditions
                n number of steps
                p steps per plotted point
        Output: t time points 
                y coordinates of roadway """
    plt.ion()  # Enable interactive plotting
    plt.figure(figsize=(8, 8))
    plt.axis([-8, 8, -8, 8])
    plt.xlabel('Roadway')
    plt.ylabel('Height')
    h = (inter[1] - inter[0]) / n
    y = np.zeros((n + 1, len(ic)))  # Initialize solution array
    y[0, :] = ic  # Set initial conditions
    t = np.zeros(n + 1)
    t[0] = inter[0]
    hw = 6  # Half-width of the roadway
    road_line, = plt.plot([], [], 'b-', linewidth=5)
    left_cable_line, = plt.plot([], [], 'r-', linewidth=1)
    right_cable_line, = plt.plot([], [], 'r-', linewidth=1)
    for k in range(n):
        for i in range(p):
            t[i + 1] = t[i] + h
            y[i + 1, :] = trapstep(f, t[i], y[i, :], h)
        y[0, :] = y[p, :]
        t[0] = t[p]
        c = hw * np.cos(y[0, 2])  # Calculate horizontal cable position
        s = hw * np.sin(y[0, 2])  # Calculate vertical cable position
        road_line.set_data([-c, c], [-s - y[0, 0], s - y[0, 0]])
        left_cable_line.set_data([-c, -c], [-s - y[0, 0], 8])
        right_cable_line.set_data([c, c], [s - y[0, 0], 8])
        plt.draw()
        plt.pause(h)
    plt.ioff()  # Disable interactive mode
    plt.show()
    return t,y

def trapstep(f, t, x, h):    # One step of the Trapezoid Method
    z1 = f(t, x)
    z2 = f(t + h, x + h*z1)
    return x + h*(z1 + z2)/2

def f(t, y):
    hw, a, W = 6, 0.2, 80
    omega = 2 * np.pi * 38/60
    a1 = np.exp(a*(y[0] - hw*np.sin(y[2])))
    a2 = np.exp(a*(y[0] + hw*np.sin(y[2])))
    z = np.zeros(4)
    z[0] = y[1]
    z[1] = -0.01*y[1]-0.4*(a1+a2-2)/a+0.2*W*np.sin(omega*t)
    z[2] = y[3]
    z[3] = -0.01*y[3]+1.2*np.cos(y[2])*(a1-a2)/(hw*a)
    return z

# Example usage
tacoma([0, 10], [1, 0, 0.001, 0], 2500, 5)
