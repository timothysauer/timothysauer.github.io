# Program 1.1 Bisection Method
def bisect(f, a, b, tol):
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