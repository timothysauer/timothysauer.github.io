import numpy as np

def newtdd(x, y):
    """
    Program 3.1 Newton Divided Differences
    Compute coefficients of interp polynomial
    Input:  x,y data points
    Output: c coefficient vector
    """
    n = len(x)  # Number of data points
    v = np.zeros((n, n))  # Create a 2D array 
    for j in range(n):
        v[j, 0] = y[j]
    for i in range(1, n):  # For column i
        for j in range(n - i):  # Fill in column 
            v[j,i] = (v[j+1,i-1]-v[j,i-1])/(x[j+i]-x[j])
    # Extract coefficients from the top of the triangle
    c = v[0, :n]
    return c

# Example usage:
# c = newtdd([0,2,3],[1,2,4])

