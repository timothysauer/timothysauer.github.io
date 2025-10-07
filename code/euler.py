import numpy as np
import matplotlib.pyplot as plt

def euler(ydot,inter, y0, n):
    """Program 6.1 Euler's Method
    Solve the IVP y' = ydot(t,y)
    Input:  ydot right-hand side function of ODE
            inter time interval of solution
            y0 initial value
            n number of time steps
    Output: t time points
            y solution"""
    t = np.zeros(n + 1)  # Time steps
    y = np.zeros(n + 1)  # Solution 
    t[0] = inter[0]
    y[0] = y0
    h = (inter[1] - inter[0])/n  # Step size
    for i in range(n):
        t[i + 1] = t[i] + h
        y[i + 1] = eulerstep(t[i], y[i], h)
    plt.plot(t, y)
    plt.xlabel('Time')
    plt.ylabel('Solution')
    plt.title('Euler Method')  
    plt.show()
    return t, y

def eulerstep(t, y, h): return y + h*ydot(t, y)

# Example usage
#def ydot(t, y): return t*y + t**3
#y = euler(ydot,[0, 1], 1., 10)
