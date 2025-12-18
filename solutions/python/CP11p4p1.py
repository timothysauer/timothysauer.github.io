import numpy as np

def mdctexample(x,n):
    le = len(x)
    M = np.zeros((n,2*n))
    for i in range(n):
        for j in range(2*n):
            M[i,j] = np.sqrt(2./n)*np.cos((i+0.5)*(j+0.5+n/2)*np.pi/n)
    N = M.T
    xout = np.array([])
    for k in range(int(le/n) - 2):
        v0 = M@x[k*n:(k+2)*n]
        v1 = M@x[(k+1)*n:(k+3)*n]
        wa = N@v0
        wb = N@v1
        w1 = wa[n:2*n]
        w2 = wb[:n]
        u1 = (w1+w2)/2
        xout = np.concatenate((xout,u1),axis = 0)
    return xout

x = np.arange(1,13)
print(mdctexample(x,4))
