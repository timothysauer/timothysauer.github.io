import numpy as np
import matplotlib.pyplot as plt

def predcorr(f, inter, ic, n, s):
    """ Program 6.8 Adams-Bashforth-Mounton second-order p-c
        Input:  f differential equation
                inter time interval
                ic initial conditions
                n number of steps
                s number of multisteps for explicit method
        Output: t time points
                y solution """
    h = (inter[1] - inter[0]) / n  # Step size
    y = np.zeros((n + 1, len(ic)))  # Initialize solution array
    fd = np.zeros((n + 1, len(ic)))  # Store derivative evaluations
    t = np.zeros(n + 1)             # Initialize time array
    y[0, :] = ic                     # Set initial condition
    t[0] = inter[0]
    fd[0,:] = f(t[0],y[0,:])
    # Start-up phase using one-step method
    for i in range(s - 1):
        t[i + 1] = t[i] + h
        y[i + 1, :] = trapstep(f, t[i], y[i, :], h)
        fd[i + 1, :] = f(t[i+1], y[i+1, :])
    # Multistep method loop
    for i in range(s - 1, n):
        t[i + 1] = t[i] + h
        y[i + 1, :] = ab2step(t[i], y[:i+1,:], fd[:i+1,:], h)  # Predictor step
        fd[i + 1, :] = f(t[i+1], y[i+1, :])
        y[i + 1, :] = am1step(t[i], y[:i+2,:], fd[:i+2,:], h)
        fd[i + 1, :] = f(t[i+1], y[i+1, :])
    return t, y

def trapstep(f, t, x, h):  # One step of the Trapezoid Method
    z1 = f(t, x)
    z2 = f(t + h, x + h*z1)
    return x + h*(z1 + z2)/2

def ab2step(t,w,fd,h):  # One step of Adams-Bashforth 2-step method 
    return w[-1,:] + h*(3*fd[-1,:]/2 - fd[-2,:]/2)

def am1step(t,w,fd,h): # One step of the AM 1-step method 
    return w[-2, :] + h*(fd[-1,:]+fd[-2,:])/2

# Example usage
def f(t, y):  return t*y + t**3
t, y = predcorr(f, [0, 1], [1], 20, 2)
plt.plot(t, y)
plt.xlabel('Time')
plt.ylabel('Solution')
plt.title('Multistep Method Solution')
plt.show()