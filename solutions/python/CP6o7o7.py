import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../code')
from exmulti import exm

def ydot(t,y): return 1 + y**2  #np.sin(y)
def yexact(t,y0): return 2.*np.arctan(np.exp(t)*np.tan(y0/2))+2*np.pi*np.floor((y0+np.pi)/(2*np.pi))
def yexact(t,y0): return np.tan(t + np.arctan(y0))
plt.rcParams.update({'font.size': 14})
inter = [0,1]
y0 = [0.5] # 100.
h = []
y4=[]
for i in range(1):  
    n = 10*2**i
    h = h+[1./n]
    t,y = exm(ydot,inter,y0,n,3)
    y4 = y4+[y[-1]]

for i in range(1):  
    n = 10*2**i
    h = h+[1./n]
    t1,y1 = exm(ydot,inter,y0,2*n,3)
    y4 = y4+[y[-1]]
ytrue = yexact(t1,y0)
y4 = np.array(y4)
#h = np.array(h)
plt.plot(t,y,t1,y1,t1,ytrue,'-k')
#plt.loglog(h,np.abs(ytrue[-1]-y4),'ks')
plt.grid(True)
#plt.xlabel('h')
#plt.ylabel('error at t = 4')
#plt.xlim(1e-3,1e-1)
#plt.ylim(1e-4,1e-2)
plt.xlim(0,1)
plt.ylim(0,10)
#plt.savefig('cp6o7o7.png')
plt.show()