def fpi(g, x0, k):
    # Initialize the array to hold values
    x = [0] * (k + 1)
    x[0] = x0
    
    for i in range(k):
        x[i + 1] = g(x[i])  # Update the next point using the function g
    
    xc = x[k]  # The approximate solution is the last computed value
    return xc

# Example usage
# Define a function g for testing. For example, g(x) = (x + 2) / 3
def g(x):
    return (x + 2) / 3

x0 = 0.0  # Starting guess
k = 10    # Number of iterations

# Compute the fixed-point iteration
solution = fpi(g, x0, k)
print("Approximate solution:", solution)
