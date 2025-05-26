import numpy as np
import matplotlib.pyplot as plt
from splinecoeff import splinecoeff

def splineplot(x, y, k):
    """Program 3.6 Cubic spline plot
    Computes and plots spline from data points
    Input: x,y arrays of data points,
       number of plotted points per segment
    Output: x1, y1 spline values at plotted points"""
    n = len(x)
    coeff = splinecoeff(x, y)
    x1 = []
    y1 = []
    for i in range(n - 1):
        # Generate equally spaced points between x[i] and x[i+1]
        xs = np.linspace(x[i], x[i + 1], k + 1)
        dx = xs - x[i]
        # Evaluate using coefficients
        ys = coeff[i, 2] * dx  # d coefficient
        ys = (ys + coeff[i, 1]) * dx  # c coefficient
        ys = (ys + coeff[i, 0]) * dx + y[i]  # b coefficient
        # Append evaluated values to the lists
        x1.extend(xs[:-1])  # Excluding the last point as it overlaps with next segment
        y1.extend(ys[:-1])  # Excluding the last point for the same reason

    # Append the last point
    x1.append(x[-1])
    y1.append(y[-1])
    
    # Plotting the original data points and the spline
    plt.plot(x, y, 'o', label='Data Points')
    plt.plot(x1, y1, label='Cubic Spline')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Cubic Spline Interpolation')
    plt.legend()
    plt.grid()
    plt.show()
    
    return np.array(x1), np.array(y1)

