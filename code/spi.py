def spi(f, r, s, t, k):
    """ Program 13.2 Successive Parabolic Interpolation
        Input:  f function
                r,s,t initial guesses
                k number of steps
        Output: x approximate minimum """
    x = [r, s, t] # Initialize points
    fr = f(r)
    fs = f(s)
    ft = f(t)
    for i in range(3, k + 3):  
        x_i=(r+s)/2-(fs-fr)*(t-r)*(t-s)/(2*((s-r)*(ft-fs)-(fs-fr)*(t-s)))
        x.append(x_i)  # Append the new point to the list
        t = s
        s = r
        r = x_i
        ft = fs
        fs = fr
        fr = f(r)  # Single function evaluation
    return x[-1]

# Example usage
def f(x): return (x - 3)**4 + 1
minimum = spi(f, 1, 1.2, 1.5, 30)
print("Approximate minimum:", minimum)
      
