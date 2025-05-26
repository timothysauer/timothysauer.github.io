import numpy as np

def romberg(f, a, b, n):
    """ Program 5.1 Romberg Integration """
    h = (b - a) / (2 ** np.arange(n))  # Step sizes
    r = np.zeros((n, n))  # Initialize Romberg tableau
    r[0, 0] = (b - a) * (f(a) + f(b)) / 2  # Initial trapezoidal rule approximation
    for j in range(1, n):
        subtotal = 0
        # Sum over the midpoint values
        for i in range(1, 2 ** (j - 1) + 1):
            subtotal += f(a + (2 * i - 1) * h[j])
        # Compute the next row of the Romberg tableau
        r[j, 0] = r[j - 1, 0] / 2 + h[j] * subtotal
        # Compute the higher order terms using Richardson Extrapolation
        for k in range(1, j + 1):
            r[j, k] = (4 ** k * r[j, k - 1] - r[j - 1, k - 1]) / (4 ** k - 1)

    return r

# Example usage
def integrand(x):
    #return np.sin(x)  # Example integrand function
    return x**3
a = 0  # Start of interval
b = np.pi  # End of interval
n = 4  # Number of rows in the Romberg tableau
a = 0
b = 1

result = romberg(integrand, a, b, n)
print("Romberg Tableau:\n", result)
	      
