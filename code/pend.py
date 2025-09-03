import numpy as np
import matplotlib.pyplot as plt

def pend(inter, ic, n):
    h = (inter[1] - inter[0]) / n  # Step size
    y = np.zeros((n + 1, 2))        # Initialize solution array
    t = np.zeros(n + 1)             # Time array
    y[0, :] = ic                     # Set initial conditions
    t[0] = inter[0]

    # Set up the plot
    plt.xlim([-1.2, 1.2])
    plt.ylim([-1.2, 1.2])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('X position')
    plt.ylabel('Y position')

    # Create line objects for bob and rod
    bob, = plt.plot([], [], 'ro', markersize=15)  # bob's position
    rod, = plt.plot([], [], 'b-', linewidth=3)   # rod's position

    # Animation loop
    for k in range(n):
        t[k + 1] = t[k] + h
        y[k + 1, :] = trapstep(t[k], y[k, :], h)  # Update state using trapezoidal step
        xbob = np.sin(y[k + 1, 0])  # x position of the bob
        ybob = -np.cos(y[k + 1, 0])  # y position of the bob
        xrod = [0, xbob]  # x coordinates of the rod
        yrod = [0, ybob]  # y coordinates of the rod
        
        # Update plot data
        rod.set_data(xrod, yrod)
        bob.set_data(xbob, ybob)
        plt.draw()
        plt.pause(h)  # Pause for h seconds

    plt.show()

def trapstep(t, x, h):
    """ One step of the Trapezoid Method """
    z1 = ydot(t, x)
    g = x + h * z1
    z2 = ydot(t + h, g)
    return x + h * (z1 + z2) / 2

def ydot(t, y):
    g = 9.81  # Acceleration due to gravity
    length = 1.0  # Length of the pendulum pivot
    z = np.zeros(2)
    z[0] = y[1]  # The first state variable is the angle theta.
    z[1] = -(g / length) * np.sin(y[0])  # The second variable is the angular velocity.
    return z

# Example usage
pend([0, 10], [np.pi / 2, 0], 200)  # Starting angle of pi/2 radians and angular velocity 0
