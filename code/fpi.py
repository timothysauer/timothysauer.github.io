def fpi(g, x0, k):
    """ Program 1.2 Fixed-Point Iteration
    Compute approximate solution of x=f(x)
    Input:  g function
            x0 initial guess
            k number of steps
    Output: xc approximate fixed point
    """
    x = [0] * (k + 1)
    x[0] = x0
    for i in range(k):
        x[i + 1] = g(x[i]) 
    xc = x[k]  
    return xc

# Example usage
# def g(x): return (x + 2) / 3
# xc = fpi(g, 0., 10)

