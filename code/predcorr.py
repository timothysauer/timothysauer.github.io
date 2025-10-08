import numpy as np
import matplotlib.pyplot as plt

def predcorr(inter, ic, n, s):
    """ Program 6.8 Adams-Bashforth-Mounton second-order p-c
        Input:  inter time interval
                ic initial conditions
                n number of steps
                s number of multisteps for explicit method
        Output: t time points
                y solution """
    h = (inter[1] - inter[0]) / n  # Step size
    y = np.zeros((n + 1, len(ic)))  # Initialize solution array
    t = np.zeros(n + 1)             # Initialize time array
    y[0, :] = ic                     # Set initial condition
    t[0] = inter[0]
    # Start-up phase using one-step method
    for i in range(s - 1):
        t[i + 1] = t[i] + h
        y[i + 1, :] = trapstep(t[i], y[i, :], h)
    # Multistep method loop
    for i in range(s - 1, n):
        t[i + 1] = t[i] + h
        f_i = ydot(t[i], y[i, :])  # Evaluate current derivative
        y[i + 1, :] = ab2step(t[i], i, y, f_i, h)  # Predictor step
        f_i_next = ydot(t[i + 1], y[i + 1, :])     # Corrector step
        y[i + 1, :] = am1step(t[i], i, y, f_i, f_i_next, h)
    # Plotting the results
    plt.plot(t, y[:, 0])
    plt.xlabel('Time')
    plt.ylabel('y')
    plt.title('Predictor-Corrector Method for ODEs')
    plt.grid()
    plt.show()
    return t, y

def trapstep(t, x, h):  # One step of the Trapezoid Method
    z1 = ydot(t, x)
    g = x + h*z1
    z2 = ydot(t + h, g)
    return x + h*(z1 + z2)/2

def ab2step(t,i,y,f,h):  # One step of AB 2-step method
    return y[i, :] + h*(3*f - ydot(t, y[i-1, :]))/2

def am1step(t,i,y,f_i,f_i_next,h): # One step of the AM 1-step method 
    return y[i, :] + h*(f_i_next + f_i)/2

def ydot(t, y):  return t * y + t**3

# Example usage
t, y = predcorr([0, 1], [1], 20, 2)
