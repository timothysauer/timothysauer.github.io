import numpy as np
import matplotlib.pyplot as plt
from splinecoeff import splinecoeff

def splineplot(x, y, k):
    """Program 3.6 Cubic spline plot
    Computes and plots spline from data points
    Input:  x,y arrays of data points,
            k number of plotted points per segment
    Output: x1, y1 spline values at plotted points"""
    n = len(x)
    coeff = splinecoeff(x, y)
    print(coeff)
    x1 = []
    y1 = []
    for i in range(n - 1):
        xs = np.linspace(x[i], x[i + 1], k + 1)
        dx = xs - x[i]
        ys = coeff[i, 2] * dx  # d coefficient
        ys = (ys + coeff[i, 1]) * dx  # c coefficient
        ys = (ys + coeff[i, 0]) * dx + y[i]  # b coefficient
        x1.extend(xs[:-1])  # Excluding the last point due to overlap
        y1.extend(ys[:-1])  
    x1.append(x[-1])
    y1.append(y[-1])
    # Plot the original data points and the spline
    plt.plot(x, y, 'o', label='Data Points')
    plt.plot(x1, y1, label='Cubic Spline')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Cubic Spline Interpolation')
    plt.legend()
    plt.grid()
    plt.show()
    return np.array(x1), np.array(y1)

# Example usage:
# x1, y1 = splineplot([1,2,3,4],[3,2,6,1],10)
