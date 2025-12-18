import numpy as np

def splinecoeffalt(x, y):
    """
    Program 3.5 Calculation of spline coefficients
    Calculates coefficients of cubic spline
    Input:  x,y vectors of data points
            plus two optional extra data vl, vr
    Output: matrix of coefficients b1,c1,d1;b2,c2,d2;..."""
    n = len(x)
    vl = 0.0  # Specify vl and vr if needed (see below)
    vr = 0.0
    A = np.zeros((n, n)) # Initialize nxn matrix A, vector r
    r = np.zeros(n)
    dx = np.diff(x)  # Define the deltas
    dy = np.diff(y)
    # Load the A matrix
    for i in range(1, n - 1):
        A[i, i - 1] = dx[i - 1]
        A[i, i] = 2 * (dx[i - 1] + dx[i])
        A[i, i + 1] = dx[i]
        r[i] = 3 * (dy[i] / dx[i] - dy[i - 1] / dx[i - 1])    
    # Set endpoint conditions (use only one of the following 5 pairs)
    #A[0, 0] = 1.  # Natural spline conditions at the first endpoint
    #A[-1, -1] = 1.  # Natural spline conditions at the last endpoint
    #A[0,0]=2;r[0]=vl;  # curvature-adj conditions (specify vl and vr)
    #A[n-1,n-1]=2;r[n-1]=vr
    #A[0,:2]=[2*dx[0],dx[0]];r[0]=3*(dy[0]/dx[0]-vl);# clamped (specify vl and vr)
    #A[n-1,n-2:]=[dx[n-2],2*dx[n-2]];r[n-1]=3*(vr-dy[n-2]/dx[n-2])
    #A[0,:2]=[1, -1]        # parabol-term conditions, for n>=3
    #A[n-1,n-2:]=[1, -1]
    A[0,:3]=[dx[1], -(dx[0]+dx[1]), dx[0]]; # not-a-knot for n>=4
    A[n-1,n-3:]=[dx[n-2], -(dx[n-3]+dx[n-2]), dx[n-3]]
    coeff = np.zeros((n, 3))  # Coefficients array
    coeff[:, 1] = np.linalg.solve(A, r)
    for i in range(n - 1):
        coeff[i,2]=(coeff[i + 1, 1] - coeff[i, 1]) / (3 * dx[i])
        coeff[i,0]=dy[i]/dx[i]-dx[i]*(2*coeff[i, 1]+coeff[i + 1, 1]) / 3
    return coeff[:-1, :3]  # Return only up to the n-1 index

# example usage:
# coeff = splinecoeff([1,2,3,4],[3,2,6,1])
