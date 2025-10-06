import numpy as np
def nest(c, x, b=None):
    """ Evaluate polynomial in nested form
    Input: c list of coefficients (constant term first)
       x argument at which to evaluate
       b list of base points, if needed
    Output: value of polynomial
    """
    d = len(c)-1
    if b is None:  # If no base points are provided
        b = np.zeros((d))
    y = c[d]  # Initialize y 
    # Iterate through coefficients in reverse order
    for i in range(d-1, -1, -1):
        y = y*(x-b[i])+c[i]  # Horner's method 
    return y
# example usage: 
# y = nest([-1,5,-3,3,2],0.5)  
