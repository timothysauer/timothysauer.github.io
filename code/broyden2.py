import numpy as np

def broyden2(f, x0, k):
    n = len(x0)
    b = np.eye(n)  # Initial approximation of the Jacobian (identity matrix)

    for i in range(k):
        # Compute the new estimate of the solution
        x = x0 - np.dot(b, f(x0))
        del_x = x - x0  # Change in x
        del_f = f(x) - f(x0)  # Change in f

        # Update Broyden's approximation of the Jacobian
        b = b + np.outer(del_x - np.dot(b, del_f), del_x) @ b / (del_x @ b @ del_f)

        x0 = x  # Update x0 for the next iteration

    return x

# Example usage
# Define a sample function, for example, f(x) = x^2 - 4
def f(x):
    return np.array([x[0]**2 - 4, x[1] - 1])  # Example function

# Initial guess and maximum steps
x0 = np.array([1.0, 2.0])
k = 10

# Call the Broyden's method function
solution = broyden2(f, x0, k)
print("Approximate solution:", solution)
