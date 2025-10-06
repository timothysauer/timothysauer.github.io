import numpy as np
import matplotlib.pyplot as plt

def newtdd(x, y):
    n = len(x)  # Number of data points
    v = np.zeros((n, n))  # Create a 2D array for the divided differences
    for j in range(n):
        v[j, 0] = y[j]
    # Fill in the divided difference table
    print(x)
    for i in range(1, n):  # For column i
        for j in range(n - i):  # Fill in column from top to bottom
            v[j, i] = (v[j + 1, i - 1] - v[j, i - 1]) / (x[j + i] - x[j])
    # Extract coefficients from the top of the triangle
    c = v[0, :n]
    return c

def nest(c, x, b=None):
    d = len(c)-1
    if b is None:  # If no base points are provided, initialize as zeros
        b = np.zeros(d)
    y = c[d]  # Initialize y with the highest degree coefficient
    # Iterate through the polynomial coefficients in reverse order
    for i in range(d-1, -1, -1):
        y = y*(x-b[i])+c[i]  # Horner's method with nested multiplication
    return y

def clickinterp():
    xl = -3
    xr = 3
    yb = -3
    yt = 3
    plt.plot([xl, xr], [0, 0], 'k', [0, 0], [yb, yt], 'k')
    plt.grid(True)
    plt.xlim(xl, xr)
    plt.ylim(yb, yt)
    xlist = []
    ylist = []
    number_of_points = 0  # Initialize counter 
    while True:
        inp = plt.ginput(1,timeout=0)  # Get one mouse click
        if len(inp) < 1:  # return key closes the program
            plt.close('all')
            break
        xnew,ynew = inp[0]
        number_of_points += 1  # k counts clicks
        xlist.append(xnew)  # Add new point to the list
        ylist.append(ynew)
        if number_of_points > 1:   # plot the interp poly if >= 2 points
            c = newtdd(xlist, ylist)  # Get interpolation coefficients
            x = np.linspace(xl, xr, 100)  # Define x coordinates of curve
            #y = newton_polynomial(c, x, xlist)  # Get y coordinates of curve
            y = nest(c,x,xlist)
            plt.clf()  # Clear the current figure
            plt.plot(xlist, ylist, 'o', x, y, 'b-')
            plt.plot([xl, xr], [0, 0], 'k', [0, 0], [yb, yt], 'k')
            plt.xlim(xl, xr)
            plt.ylim(yb, yt)
            plt.grid(True)
            plt.draw()
    plt.show()

# Run the interactive interpolation
clickinterp()
