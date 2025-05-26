import numpy as np
import numpy as np
from nest import nest
from newtdd import newtdd

def sin1(x):
    """ 
    Program 3.3 Building a sin calculator key, attempt #1
    Approximates sin curve with degree 3 polynomial
    (Caution: do not use to build bridges,
    at least until we have discussed accuracy.)
    Input: x
    Output: approximation for sin(x) """
    b = np.pi * np.array([0, 1/6, 2/6, 3/6])    # Base points
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

