import numpy as np
def nest(c, x, b=None):
    d = len(c)-1
    if b is None:  # If no base points are provided, initialize as zeros
        b = np.zeros(d)
    y = c[d]  # Initialize y with the highest degree coefficient
    # Iterate through the polynomial coefficients in reverse order
    for i in range(d-1, -1, -1):
        y = y*(x-b[i])+c[i]  # Horner's method with nested multiplication
    return y
      
