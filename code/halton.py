import numpy as np

def halton(p, n):
    """ Program 9.1 Quasi-random number generator
        Input:  p prime number
                n requested number 
        Output: u array of n quasi-random numbers """
    # Find the largest number of digits in the base p representation
    b = np.zeros(int(np.ceil(np.log(n)/np.log(p))), dtype=int)
    u = np.zeros(n)  # Initialize array for quasi-random numbers
    for j in range(n):
        i = 0
        b[0] += 1  # increment the current integer
        while b[i] > p - 1:  # Loop for carrying in base p
            b[i] = 0
            i += 1
            b[i] += 1
        u[j] = 0
        for k in range(len(b)):  # Sum up the reversed digits, and use the ...
            u[j] += b[k]*p**(-k-1)  # ... fractional representation
    return u

# Example usage
u = halton(2, 100)

