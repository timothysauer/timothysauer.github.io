import numpy as np

def gss(f, a, b, k):
    """ Program 13.1 Golden Section Seach for minimum
        Start with unimodal f and minimum in [a,b]
        Input:  f unimodal function to minimize
                a, b endpoints of interval
                k number of steps
        Output: x approximate minimum """
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
    return (a + b)/2  # Approximate minimum location

# Example usage
def f(x): return (x - 3)**2
minimum = gss(f, 0, 5, 20)
print("Approximate minimum:", minimum)


