def halton(p, n):
    """ Program 9.1 Quasi-random number generator
    Input:  p prime number
            n requested quantity of numbers 
    Output: u array of n quasi-random numbers """
    u = []  # Initialize list for quasi-random numbers
    for k in range(1,n+1):
        j = k
        num = 0
        i = 1
        while j > 0:
            num += (j % p)/p**i
            j //= p
            i += 1
        u.append(num)
    return u

# Example usage
u = halton(2, 100)


