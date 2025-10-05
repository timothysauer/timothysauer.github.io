import numpy as np

def gss(f, a, b, k):
    g = (np.sqrt(5) - 1)/2
    x1 = a + (1 - g)*(b - a)
    x2 = a + g*(b - a)
    f1 = f(x1)
    f2 = f(x2)

    for _ in range(k):
        if f1 < f2:  # If f(x1) < f(x2), replace b with x2
            b = x2
            x2 = x1
            x1 = a + (1 - g)*(b - a)
            f2 = f1
            f1 = f(x1)  # Single function evaluation
        else:  # Otherwise, replace a with x1
            a = x1
            x1 = x2
            x2 = a + g*(b - a)
            f1 = f2
            f2 = f(x2)  # Single function evaluation

    y = (a + b)/2  # Approximate minimum location
    return y

# Example usage
# Define a sample function, for example f(x) = (x - 2)**2
def f(x):
    return (x - 2)**2

# Specify the interval and number of steps
a = 0  # Left endpoint
b = 5  # Right endpoint
k = 20  # Number of iterations

# Compute the minimum
minimum = gss(f, a, b, k)
print("Approximate minimum:", minimum)
