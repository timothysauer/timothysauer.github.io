import numpy as np

def romberg(f, a, b, n):
    """ Program 5.1 Romberg Integration 
        Compute approximation to definite integral
        Input:  f function (integrand)
                a,b integration interval
                n number of rows of Romberg tableau
        Output: r Romberg tableau """
    h = (b - a)/(2 ** np.arange(n))  # Step sizes
    r = np.zeros((n, n))  # Initialize Romberg tableau
    r[0, 0] = (b - a)*(f(a) + f(b))/2  # Trapezoidal rule
    for j in range(1, n):
        subtotal = 0
        for i in range(1, 2 ** (j - 1) + 1):
            subtotal += f(a + (2 * i - 1) * h[j])
        r[j, 0] = r[j - 1, 0] / 2 + h[j] * subtotal
        for k in range(1, j + 1):    #  Richardson Extrapolation
            r[j, k] = (4**k*r[j, k-1] - r[j-1, k-1])/(4**k-1)
    return r

# Example usage
# def integrand(x): return np.log(x)
# result = romberg(integrand, 1, 2, 4)

	      
