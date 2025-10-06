import numpy as np

def broyden2(f, x0, k):
    """Program 2.3 Broyden's Method II
    Solve nonlinear system of equations
    Input:  x0 initial vector
            k maximum number of steps
    Output: x approximate solution
    """
    n = len(x0)
    b = np.eye(n)  # Initial approximation of the Jacobian 
    for i in range(k):
        x = x0 - np.dot(b, f(x0)) # New solution estimate
        del_x = x - x0  # Change in x
        del_f = f(x) - f(x0)  # Change in f
        b = b+np.outer(del_x-np.dot(b, del_f),del_x)@b/(del_x@b@del_f)
        x0 = x  # Update x0 for the next iteration
    return x

# Example usage
# Define a sample function, for example, f(x) = x^2 - 4
# def f(x): return np.array([x[0]**2 - 4, x[1] - 1])  
# x = broyden2(f, np.array([1.0, 2.0]), 10)
