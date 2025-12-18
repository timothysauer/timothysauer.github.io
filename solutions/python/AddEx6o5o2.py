import numpy as np
from scipy.integrate import solve_ivp
import math
import matplotlib.pyplot as plt
import sys
sys.path.append('../../code')

def ydot(t,y): return (1. + y**2)/2.
def yexact(t,y0): return 1./np.cos(t+np.pi/6)+np.tan(t+np.pi/6)
plt.rcParams.update({'font.size': 14})
inter = [0,1]
y0 = np.sqrt(3) # 100.
h = []
plt.figure()
ax = plt.axes()
#t,y = trap(ydot,inter,y0,100)
#ytrue = yexact(t,y0)
#plt.loglog(h,np.abs(ytrue[-1]-y4),'ks')
#ax.plot(t,ytrue,linewidth=3)

#y = solve_ivp(ydot,[0,1],[y0],method='RK45',rtol=1e-8)
y = solve_ivp(ydot,[0,1],[y0],method='BDF',rtol=1e-8)
print(len(y.t))
print(y.t,y.y)
ax.plot(y.t,np.squeeze(y.y),'r-')

ax.grid(True)
ax.set_xlabel('t')
ax.set_ylabel('y')
#plt.xlim(1e-3,1e-1)
#plt.ylim(1e-4,1e-2)
#plt.savefig('cp6o2o9a.png')
plt.show()