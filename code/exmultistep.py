import numpy as np
import matplotlib.pyplot as plt

def exmultistep(f, inter, ic, n, s):
    """ Program 6.7 Multistep Method
        Input:  inter time interval
                ic initial conditions
                n number of steps
                s number of multisteps, e.g. 2 for 2-step method
        Output: t time points
                y solution of IVP """
    h = (inter[1] - inter[0]) / n
    y = np.zeros((n + 1, len(ic)))  # Initialize solution array
    fd = np.zeros((n + 1, len(ic)))  # Store derivative evaluations
    t = np.zeros(n + 1)  # Initialize time steps
    y[0, :] = ic     
    t[0] = inter[0]
    fd[0,:] = f(t[0],y[0,:])
    for i in range(s - 1):  # Start-up phase, using one-step method
        t[i+1] = t[i] + h
        y[i+1, :] = trapstep(f, t[i], y[i, :], h)
        fd[i+1, :] = f(t[i+1], y[i+1, :])
    for i in range(s - 1, n):  # multistep method loop
        t[i+1] = t[i] + h
        y[i+1, :] = ab2step(t[i],y[:i+1,:],fd[:i+1,:],h)
        fd[i+1, :] = f(t[i+1], y[i+1, :])
    return t, y

def trapstep(f, t, x, h):  # One step of the Trapezoid Method
    z1 = f(t, x)  
    z2 = f(t + h, x + h*z1)
    return x + h*(z1 + z2)/2

def ab2step(t,w,fd,h):  # One step of Adams-Bashforth 2-step method 
    return w[-1,:] + h*(3*fd[-1,:]/2 - fd[-2,:]/2)

def unstable2step(t,w,fd,h):  # One step of an unstable 2-step method
    return -w[-1, :] + 2*w[-2,:] + h*(5*fd[-1, :]/2 +fd[-2, :]/2)

def weaklystable2step(t,w,fd,h): # One step of weakly-stable 2-step method
    return w[-2, :] + h*2*fd[-1, :]

# Example usage
def f(t, y):  return t*y + t**3
t, y = exmultistep(f, [0, 1], [1], 20, 2)
plt.plot(t, y)
plt.xlabel('Time')
plt.ylabel('Solution')
plt.title('Multistep Method Solution')
plt.show()