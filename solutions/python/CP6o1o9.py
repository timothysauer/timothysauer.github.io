import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../code')
from euler import euler

def ydot(t,y): return np.sin(y)
def yexact(t,y0): return 2.*np.arctan(np.exp(t)*np.tan(y0/2))+2*np.pi*np.floor((y0+np.pi)/(2*np.pi))
plt.rcParams.update({'font.size': 14})
inter = [0,4]
y0 = 100. # 100.
h = []
y4=[]
for i in range(6):  
    n = 40*2**i
    h = h+[4./n]
    t,y = euler(ydot,inter,y0,n)
    y4 = y4+[y[-1]]
ytrue = yexact(t,y0)
y4 = np.array(y4)
h = np.array(h)
plt.loglog(h,np.abs(ytrue[-1]-y4),'ks')
plt.grid(True)
plt.xlabel('h')
plt.ylabel('error at t = 4')
plt.xlim(1e-3,1e-1)
plt.ylim(1e-4,1e-2)
plt.savefig('cp6_1_9d.png')
plt.show()