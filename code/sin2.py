import numpy as np
from nest import nest
from newtdd import newtdd

def sin2(x):
    """ 
    Program 3.4 Building a sin calculator key, attempt #2
    Approximates sin curve with degree n-1 polynomial
    Input:  x
    Output: y approximation for sin(x) """
    n = 10  # Degree of the polynomial
    b = np.pi/4 + (np.pi/4)*np.cos((np.arange(1,2*n,2)*np.pi)/(2*n))  # Chebyshev
    yb = np.sin(b)  # Corresponding sine values
    c = newtdd(b, yb)  # Get coefficients using Newton divided difference
    # Adjust x to the fundamental domain
    s = 1  # Correct the sign of sin
    x1 = x % (2 * np.pi)  # Normalize x to [0, 2*pi]
    if x1 > np.pi:
        x1 = 2 * np.pi - x1
        s = -1
    if x1 > np.pi / 2:
        x1 = np.pi - x1
    y = s * nest(c, x1, b) # Evaluate the interpolating polynomial
    return y

# Example usage:
# y = sin2(np.pi/4)
