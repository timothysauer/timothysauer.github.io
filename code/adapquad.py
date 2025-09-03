import numpy as np

def adapquad(f, a0, b0, tol0):
    sum = 0
    n = 1
    a = [a0]
    b = [b0]
    tol = [tol0]
    app = [trap(f, a[0], b[0])]

    while n > 0:  # n is the current position at the end of the list
        c = (a[n-1] + b[n-1]) / 2  # Adjusted to 0-based index
        oldapp = app[n-1]
        
        app[n-1] = trap(f, a[n-1], c)
        app.append(trap(f, c, b[n-1]))  # Append new value for the second half

        if abs(oldapp - (app[n-1] + app[n])) < 3 * tol[n-1]:
            sum += app[n-1] + app[n]  # Success
            n -= 1  # Done with this interval
        else:  # Divide into two intervals
            b.append(b[n-1])  # Keep ending point
            b[n-1] = c  # New end for the first interval
            a.append(c)  # New start for the second interval
            tol[n-1] = tol[n-1] / 2  # Reduce tolerance
            tol.append(tol[n-1])  # Same tolerance for the new interval
            n += 1  # Go to end of list, repeat

    return sum

def trap(f, a, b):
    return (f(a) + f(b)) * (b - a) / 2

# Example usage
# Define a function to integrate, for example f(x) = x**2
def f(x):
    return x ** 2

# Call adapquad function
integral = adapquad(f, 0, 1, 0.001)
print("Approximate integral:", integral)