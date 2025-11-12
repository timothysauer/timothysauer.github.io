import numpy as np
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

def ydot(t,y): return np.sin(y)
def yexact(t,y0): return 2.*np.arctan(np.exp(t)*np.tan(y0/2))+2*np.pi*np.floor((y0+np.pi)/(2*np.pi))
plt.rcParams.update({'font.size': 14})
inter = [0,4]
y0 = 100. 
h = []
y4=[]
for i in range(1):  
    n = 40*2**i
    h = h+[4./n]
    t,y = trap(ydot,inter,y0,n)
    y4 = y4+[y[-1]]
ytrue = yexact(t,y0)
y4 = np.array(y4)
h = np.array(h)
plt.plot(t,ytrue,t,y)
#plt.loglog(h,np.abs(ytrue[-1]-y4),'ks')
plt.grid(True)
plt.xlabel('t')
plt.ylabel('y')
plt.savefig('cp6o2o7b.png')
plt.show()