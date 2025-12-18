import numpy as np
import matplotlib.pyplot as plt

def predcorr(ydot, inter, ic, n, s):
    """ Program 6.8 Adams-Bashforth-Mounton second-order p-c
        Input:  ydot differential equation
                inter time interval
                ic initial conditions
                n number of steps
                s number of multisteps for explicit method
        Output: t time points
                y solution """
    h = (inter[1] - inter[0]) / n  # Step size
    y = np.zeros((n + 1, len(ic)))  # Initialize solution array
    f = np.zeros((n + 1, len(ic)))  # Store derivative evaluations
    t = np.zeros(n + 1)             # Initialize time array
    y[0, :] = ic                     # Set initial condition
    t[0] = inter[0]
    f[0,:] = ydot(t[0],y[0,:])
    # Start-up phase using one-step method
    for i in range(s - 1):
        t[i + 1] = t[i] + h
        y[i + 1, :] = trapstep(t[i], y[i, :], h)
        f[i + 1, :] = ydot(t[i+1], y[i+1, :])
    # Multistep method loop
    for i in range(s - 1, n):
        t[i + 1] = t[i] + h
        y[i + 1, :] = ab2step(t[i], y[:i+1,:], f[:i+1,:], h)  # Predictor step
        f[i + 1, :] = ydot(t[i+1], y[i+1, :])
        y[i + 1, :] = am1step(t[i], y[:i+2,:], f[:i+2,:], h)
        f[i + 1, :] = ydot(t[i+1], y[i+1, :])
    return t, y

def trapstep(ydot, t, x, h):  # One step of the Trapezoid Method
    z1 = ydot(t, x)
    g = x + h*z1
    z2 = ydot(t + h, g)
    return x + h*(z1 + z2)/2

def ab2step(t,w,f,h):  # One step of Adams-Bashforth 2-step method 
    return w[-1,:] + h*(3*f[-1,:]/2 - f[-2,:]/2)

def am1step(t,w,f,h): # One step of the AM 1-step method 
    return w[-2, :] + h*(f[-1,:]+f[-2,:])/2

# Example usage
#def ydot(t, y):  return t*y + t**3
#t, y = predcorr(ydot, [0, 1], [1], 20, 2)
#plt.plot(t, y)
#plt.xlabel('Time')
#plt.ylabel('Solution')
#plt.title('Multistep Method Solution')
#plt.show()