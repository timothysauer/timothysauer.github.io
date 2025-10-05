def spi(f, r, s, t, k):
    # Initialize points
    x = [r, s, t]
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
    return x

# Example usage
def example_function(x):
    return (x - 2)**4 + 1  # Parabola with vertex at (2, 1)
initial_guesses = (1, 1.2, 1.5)  # Initial guesses r, s, t
steps = 20  # Number of iterations
result = spi(example_function, *initial_guesses, steps)
print("Approximated minima (x):", result)
      
