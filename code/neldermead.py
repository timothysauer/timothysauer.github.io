import numpy as np

def nelder_mead(f, xbar, rad, k):
    n = len(xbar)
    # Initialize the simplex vertices
    x = np.zeros((n, n + 1))
    x[:, 0] = xbar
    x[:, 1:n + 1] = xbar[:, np.newaxis] + rad * np.eye(n)
    y = np.zeros(n + 1) # Evaluate at each vertex of the simplex
    for j in range(n + 1):
        y[j] = f(x[:, j])
    y, r = zip(*sorted(zip(y, range(n + 1)))) # Sort function values 
    y = np.array(y)                 # and rearrange simplex vertices
    x = x[:, list(r)]
    for i in range(k):
        xbar = np.mean(x[:, :n], axis=1)  # Centroid of the simplex 
        xh = x[:, n]  # Worst vertex
        xr = 2 * xbar - xh   # Reflection
        yr = f(xr)
        if yr < y[n]:  # If reflection is better than the worst point
            if yr < y[0]:  # If reflection is better than the best point
                # Expansion
                xe = 3 * xbar - 2 * xh
                ye = f(xe)
                if ye < yr:  # Accept expansion
                    x[:, n] = xe
                    y[n] = ye
                else:  # Accept reflection
                    x[:, n] = xr
                    y[n] = yr
            else:  # Accept reflection
                x[:, n] = xr
                y[n] = yr
        else:  # Reflection is not better than the worst point
            if yr < y[n]:  # Try outside contraction
                xoc = 1.5 * xbar - 0.5 * xh
                yoc = f(xoc)
                if yoc < yr:  # Accept outside contraction
                    x[:, n] = xoc
                    y[n] = yoc
                else:  # Shrink simplex toward the best point
                    for j in range(1, n + 1):
                        x[:, j] = 0.5 * x[:, 0] + 0.5 * x[:, j]
                        y[j] = f(x[:, j])
            else:  # Try inside contraction
                xic = 0.5 * xbar + 0.5 * xh
                yic = f(xic)
                if yic < y[n]:  # Accept inside contraction
                    x[:, n] = xic
                    y[n] = yic
                else:  # Shrink simplex toward the best point
                    for j in range(1, n + 1):
                        x[:, j] = 0.5 * x[:, 0] + 0.5 * x[:, j]
                        y[j] = f(x[:, j])
        
        # Resort the function values and rank the vertices
        y, r = zip(*sorted(zip(y, range(n + 1))))
        y = np.array(y)
        x = x[:, list(r)]
    return x, y

# Example usage
# Define a simple test function, for example:
f = lambda x:5*x[0]**4+4*x[0]**2*x[1]-x[0]*x[1]**3+4*x[1]**4-x[0]
x_init = np.array([0, 0])  # Initial guess
radius = 1  # Initial search radius
steps = 30  # Number of iterations
x_final, y_final = nelder_mead(f, x_init, radius, steps)
print("Final vertices of the simplex:\n", x_final)
print("Function values:\n", y_final)
