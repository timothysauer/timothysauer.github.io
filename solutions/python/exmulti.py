import numpy as np
import matplotlib.pyplot as plt


def exm(ydot, inter, ic, n, s):
    h = (inter[1] - inter[0]) / n
    y = np.zeros((n + 1, len(ic)))  # Initialize solution array
    f = np.zeros((n + 1, len(ic)))  # Store derivative evaluations
    t = np.zeros(n + 1)  # Initialize time steps
    y[0, :] = ic     
    t[0] = inter[0]
    f[0,:] = ydot(t[0],y[0,:])
    for i in range(s - 1):  # Start-up phase, using one-step method
        t[i+1] = t[i] + h
        y[i+1, :] = rk4step(ydot, t[i], y[i, :], h)
        f[i+1, :] = ydot(t[i+1], y[i+1, :])
    for i in range(s - 1, n):  # multistep method loop
        t[i+1] = t[i] + h
        y[i+1, :] = ab2step(t[i],y[:i+1,:],f[:i+1,:],h)
        f[i+1, :] = ydot(t[i+1], y[i+1, :])
    return t, y

def trapstep(ydot, t,x,h):
    z1 = ydot(t, x)  # One step of the Trapezoid Method
    g = x + h*z1
    z2 = ydot(t + h, g)
    return x + h*(z1 + z2)/2

def rk4step(ydot, t, w, h):
    """ One step of the Runge-Kutta order 4 method """
    s1 = ydot(t, w)
    s2 = ydot(t + h / 2., w + h * s1 / 2.)
    s3 = ydot(t + h / 2., w + h * s2 / 2.)
    s4 = ydot(t + h, w + h * s3)
    return w + h * (s1 + 2. * s2 + 2. * s3 + s4) / 6.

def ab2step(t,w,f,h):  # One step of Adams-Bashforth 2-step method 
    return w[-1,:] + h*(3*f[-1,:]/2 - f[-2,:]/2)

def ab3step(t,w,f,h):  # One step of Adams-Bashforth 2-step method 
    return w[-1,:] + h*(23*f[-1,:] - 16*f[-2,:]+5*f[-3,:])/12


# Example usage
#t, y = exmultistep([0, 1], [1], 20, 2)