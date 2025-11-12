import numpy as np
import math
import matplotlib.pyplot as plt
import sys
sys.path.append('../../code')

def trap(ydot,inter, y0, n):
    t = np.zeros(n + 1)  # Time steps
    y = np.zeros(n + 1)  # Solution 
    t[0] = inter[0]
    y[0] = y0
    h = (inter[1] - inter[0])/n  # Step size
    for i in range(n):
        t[i + 1] = t[i] + h
        y[i + 1] = trapstep(ydot, t[i], y[i], h)
    return t,y
def trapstep(ydot, t, y, h): return y + h/2.*(ydot(t, y)+ydot(t+h,y+h*ydot(t,y)))
def ydot(t,y): return (1. + y**2)/2.
def yexact(t,y0): return 1./np.cos(t+np.pi/6*0)+np.tan(t+np.pi/6*0)
plt.rcParams.update({'font.size': 14})
inter = [0,1]
y0 = 1. #np.sqrt(3) # 100.
h = []
plt.figure()
ax = plt.axes()
t,y = trap(ydot,inter,y0,100)
ytrue = yexact(t,y0)
#plt.loglog(h,np.abs(ytrue[-1]-y4),'ks')
ax.plot(t,ytrue,linewidth=3)
for i in range(6):
    n = 10*2**i
    h = h+[1./n]
    t,y = trap(ydot,inter,y0,n)
    ax.plot(t,y,'r-')

ax.grid(True)
ax.set_xlabel('t')
ax.set_ylabel('y')
#plt.xlim(1e-3,1e-1)
#plt.ylim(1e-4,1e-2)
plt.savefig('cp6o2o9a.png')
plt.show()