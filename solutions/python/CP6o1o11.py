import numpy as np
import math
import matplotlib.pyplot as plt
import sys
sys.path.append('../../code')
from euler import euler

def ydot(t,y): return (1. + y**2)/2.
def yexact(t,y0): return 1./np.cos(t+np.pi/6*0)+np.tan(t+np.pi/6*0)
plt.rcParams.update({'font.size': 14})
inter = [0,1]
y0 = 1.#np.sqrt(3) # 100.
h = []
plt.figure()
ax = plt.axes()
t,y = euler(ydot,inter,y0,100)
ytrue = yexact(t,y0)
#plt.loglog(h,np.abs(ytrue[-1]-y4),'ks')
ax.plot(t,ytrue,linewidth=3)
for i in range(6):
    n = 10*2**i
    h = h+[1./n]
    t,y = euler(ydot,inter,y0,n)
    ax.plot(t,y,'r-')

ax.grid(True)
ax.set_xlabel('t')
ax.set_ylabel('y')
#plt.xlim(1e-3,1e-1)
#plt.ylim(1e-4,1e-2)
plt.savefig('cp6_1_11a.png')
plt.show()