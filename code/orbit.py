import numpy as np
import matplotlib.pyplot as plt

def orbit(inter, ic, n, p):
    """
    Input:  inter time interval
            ic initial conditions [x0,vx0,y0,vy0]
            n number of steps
            p step per plotted point
    Output: t time points
            traj orbital trajectory
    """
    h = (inter[1] - inter[0]) / n  # Step size
    x0, vx0, y0, vy0 = ic  # Initial conditions
    traj = np.zeros((n+1, 4))  # Initialize position and velocity array
    y = np.zeros((p+1, 4))
    y[0, :] = [x0, vx0, y0, vy0]
    t = np.zeros((n+1))  # Time array
    t[0] = inter[0]
    plt.figure()
    plt.xlim(-5, 5)
    plt.ylim(-5, 5)
    plt.plot(0, 0, 'yo', markersize=16) # Draw sun at origin
    # Prepare line objects for head and tail of the orbit
    head, = plt.plot([], [], 'ro', markersize=8)
    tail, = plt.plot([], [], 'b-', linestyle='-', markersize=4)
    plt.ion()  # Enable interactive mode for real-time plotting
    for k in range(n // p):
        for i in range(p):
            t[k*p+i+1] = t[k*p+i] + h
            y[i + 1, :] = eulerstep(t[i], y[i, :], h)
        traj[k*p+1:(k+1)*p+1,:] = y[1:,:]
        # Update the position for the next frame
        set_position(head, y[p, 0], y[p, 2])  # Update head position
        set_position(tail, traj[1:(k+1)*p+1, 0], traj[1:(k+1)*p+1, 2]) 
        y[0, :] = y[p, :]  # Update last computed position
        plt.draw()
        plt.pause(0.03)  # Pause to visualize updates
    plt.ioff()  # Disable interactive mode
    plt.show()  # Show the plot window
    return t, traj

def eulerstep(t, x, h):
    return x + h * ydot(t, x)

def ydot(t, x):
    m2 = 3  # Mass of the second object
    g = 1   # Gravitational acceleration
    mg2 = m2*g
    px1, py1, vx1, vy1 = x[0], x[2], x[1], x[3]  # Unpack values
    px2, py2 = 0, 0  # Position of the sun located at the origin
    dist = np.sqrt((px2 - px1)**2 + (py2 - py1)**2)
    z = np.zeros(4)
    z[0] = vx1
    z[1] = (mg2*(px2 - px1))/(dist**3)  # Acceleration in x
    z[2] = vy1
    z[3] = (mg2*(py2 - py1))/(dist**3)  # Acceleration in y
    return z

def set_position(line, xdata, ydata):
    line.set_data(xdata, ydata)  # Function to set the position of plot points

# Example usage
orbit([0, 100], [0, 1, 2, 0], 10000, 5)
