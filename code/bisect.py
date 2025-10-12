def bisect(f, a, b, tol):
    """ Program 1.1 Bisection Method
    Compute approximate solution of f(x)=0
    Input:  f function
            a,b endpoints of interval
            tol error tolerance
    Output: xc approximate solution
    """
    # Check if the conditions for bisection method are met
    if f(a) * f(b) >= 0:
        raise ValueError('f(a) * f(b) must be less than 0!')
    fa = f(a)
    fb = f(b)
    while (b - a) / 2 > tol:
        c = (a + b) / 2
        fc = f(c)
        if fc == 0:  # c is a solution, done
            break
        if (fa * fc) < 0:  # a and c make the new interval
            b = c
            fb = fc
        else:  # c and b make the new interval
            a = c
            fa = fc
    xc = (a + b) / 2  # New midpoint is the best estimate
    return xc

# Example usage:
def f(x): return x**3 + x - 1
xc = bisect(f, 0, 1, 0.5e-4)