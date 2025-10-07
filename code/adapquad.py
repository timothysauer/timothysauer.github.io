import numpy as np

def adapquad(f, a0, b0, tol0):
    """Program 5.2 Adaptive Quadrature
    Compute approximation to definite integral
    Input:  f function (integrand)
            a0, b0 integration interval
            TOL error tolerance
    Output: sum approximate integral"""
    s = 0
    n = 1
    a, b = [a0],[b0]
    tol = [tol0]
    app = [trap(f, a[0], b[0])]
    while n > 0:  # n is the current position at the end of the list
        c = (a[n-1] + b[n-1])/2 
        oldapp = app[n-1]
        app[n-1] = trap(f, a[n-1], c)
        app.append(trap(f, c, b[n-1]))  # Second half
        if abs(oldapp - (app[n-1] + app[n])) < 3*tol[n-1]:
            s += app[n-1] + app[n]  # Success
            n -= 1  # Done with this interval
            a.pop() # Remove interval [a,b] from list
            b.pop()
            app.pop() # Remove approximations of both half-intervals
            app.pop()
            tol.pop() # Remove tolerance of completed interval
        else:  # Divide into two intervals
            a.append(c)  # Add the new interval [c,b] to list
            b.append(b[n-1])  
            b[n-1] = c  # Replace [a,b] with [a,c] 
            tol[n-1] = tol[n-1]/2  # Reduce tolerance
            tol.append(tol[n-1])  # Same tolerance for the new interval
            n += 1  # Go to end of list, repeat
    return s

def trap(f, a, b):
    return (f(a) + f(b))*(b - a)/2

# Example usage
#def f(x): return np.sin(x) 
#integral = adapquad(f, 0, np.pi, 0.0001)