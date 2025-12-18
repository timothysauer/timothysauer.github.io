import numpy as np
import matplotlib.pyplot as plt
from nest import nest
from newtdd import newtdd

def clickinterp():
    """
    Program 3.2 Interactive Polynomial Interpolation 
    Plot interpolating polynomial thru input points
    Input:  x,y from mouse clicks
    Output: plot of interpolating polynomial
    """
    xl,xr,yb,yt  = -3, 3, -3, 3
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
