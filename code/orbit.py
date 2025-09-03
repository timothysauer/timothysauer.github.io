import numpy as np
import matplotlib.pyplot as plt

def orbit(inter, ic, n, p):
    h = (inter[1] - inter[0]) / n  # Step size
    x0, vx0, y0, vy0 = ic  # Initial conditions
    y = np.zeros((n // p + 1, 4))  # Initialize position and velocity array
    y[0, :] = [x0, vx0, y0, vy0]
    t = np.zeros(n // p + 1)  # Time array
    t[0] = inter[0]

    # Set up the plot
    plt.figure()
    plt.xlim(-5, 5)
    plt.ylim(-5, 5)
    plt.xlabel('X Position')
    plt.ylabel('Y Position')

    # Draw sun as a yellow point at the origin
    plt.plot(0, 0, 'yo', markersize=15)

    # Prepare line objects for head and tail of the orbit
    head, = plt.plot([], [], 'ro', markersize=5)
    tail, = plt.plot([], [], 'b-', linestyle='-', markersize=2)

    plt.ion()  # Enable interactive mode for real-time plotting

    for k in range(n // p):
        for i in range(p):
            t[i + 1] = t[i] + h
            y[i + 1, :] = eulerstep(t[i], y[i, :], h)
        
        # Update the position for the next frame
        y[0, :] = y[p, :]  # Update the first position to the last computed position
        set_position(head, y[0, 0], y[0, 2])  # Update head position
        set_position(tail, y[1:p, 0], y[1:p, 2])  # Update tail positions
        
        plt.draw()
        plt.pause(0.01)  # Pause to visualize updates

    plt.ioff()  # Disable interactive mode
    plt.show()  # Show the plot window

def eulerstep(t, x, h):
    """One step of the Euler method"""
    return x + h * ydot(t, x)

def ydot(t, x):
    m2 = 3  # Mass of the second object
    g = 1   # Gravitational acceleration
    mg2 = m2 * g
    px1, py1, vx1, vy1 = x[0], x[2], x[1], x[3]  # Unpack values
    px2, py2 = 0, 0  # Position of the sun located at the origin
    dist = np.sqrt((px2 - px1) ** 2 + (py2 - py1) ** 2)

    # Derivatives
    z = np.zeros(4)
    z[0] = vx1
    z[1] = (mg2 * (px2 - px1)) / (dist ** 3)  # Acceleration in x
    z[2] = vy1
    z[3] = (mg2 * (py2 - py1)) / (dist ** 3)  # Acceleration in y

    return z

def set_position(line, xdata, ydata):
    """Function to set the position of plot points."""
    line.set_data(xdata, ydata)

# Example usage
orbit([0, 100], [0, 1, 2, 0], 10000, 5)
