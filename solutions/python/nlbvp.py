import numpy as np
import matplotlib.pyplot as plt

def nlbvpfd(inter, bv, n):
    """ Program 7.1 Nonlinear Finite Difference for BVP
        Input:  inter time interval
                bv boundary values
                n number of steps
        Output: w solution """
    a, b = inter[0], inter[1]
    ya, yb = bv[0], bv[1]
    h = (b - a)/(n + 1)  # Step size
    w = np.zeros(n)  # Initialize solution array w
    for i in range(20):  # Loop for Newton's method
        w = w - np.linalg.solve(jac(w,inter,bv,n),f(w,inter,bv,n))
    t = [a]+list(a+(1+np.arange(n))*h)+[b]
    w  = [ya]+list(w)+[yb]
    return t, w

def f(w, inter, bv, n):
    y = np.zeros(n)
    h = (inter[1] - inter[0])/(n + 1)
    y[0] = bv[0] - (2 + h**2)*w[0] + h**2*w[0]**2 + w[1]
    y[-1] = w[-2] - (2 + h**2)*w[-1] + h**2*w[-1]**2 + bv[1]
    for i in range(1, n-1):
        y[i] = w[i-1] - (2+h**2)*w[i]+h**2*w[i]**2+w[i+1]
    return y

def jac(w, inter, bv, n):
    a = np.zeros((n, n))
    h = (inter[1] - inter[0])/(n + 1)
    for i in range(n):
        a[i, i] = 2*h**2*w[i] - 2 - h**2
    for i in range(n - 1):
        a[i, i + 1] = 1
        a[i + 1, i] = 1
    return a

# Example usage
t, y = nlbvpfd([0, 1], [1, 4], 40)
plt.plot(t,y)
plt.xlabel('t')
plt.ylabel('y')
plt.title('Nonlinear BVP Solution')
plt.grid(True)
plt.show()