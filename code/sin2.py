import numpy as np

def newtdd(x, y):
    """ Compute coefficients of the Newton divided difference polynomial. """
    n = len(x)
    v = np.zeros((n, n))
    
    # Fill the first column with y values
    v[:, 0] = y
    
    # Fill in the divided difference table
    for j in range(1, n):
        for i in range(n - j):
            v[i, j] = (v[i + 1, j - 1] - v[i, j - 1]) / (x[i + j] - x[i])
    
    # The coefficients are in the first row
    return v[0, :]

def nest(d, c, x, b):
    """ Evaluate polynomial at x using nested multiplication (Horner's method). """
    y = c[d]  # Start with the highest degree coefficient
    for i in range(d - 1, -1, -1):
        y = y * (x - b[i]) + c[i]
    return y

import numpy as np
from nest import nest
from newtdd import newtdd

def sin2(x):
    """ 
    Program 3.4 Building a sin calculator key, attempt #2
    Approximates sin curve with degree n-1 polynomial
    Input: x
    Output: approximation for sin(x) """
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
