import numpy as np

def newtdd(x, y):
    n = len(x)  # Number of data points
    v = np.zeros((n, n))  # Create a 2D array for the divided differences
    for j in range(n):
        v[j, 0] = y[j]
    # Fill in the divided difference table
    for i in range(1, n):  # For column i
        for j in range(n - i):  # Fill in column from top to bottom
            v[j, i] = (v[j + 1, i - 1] - v[j, i - 1]) / (x[j + i] - x[j])
    # Extract coefficients from the top of the triangle
    c = v[0, :n]
    return c

