import numpy as np
import matplotlib.pyplot as plt
from nest import nest
from newtdd import newtdd

y = np.array([91.69,92.06,92.56,95.03,95.11,88.94,90.25,94.31,96.33,100.04])
t = np.array([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])+2015.
t1 = np.linspace(2015.,2030.,400)
A = np.array([np.ones((10)),t]).T
c = np.linalg.solve(A.T@A,A.T@y)
RMSE = np.linalg.norm(y-A@c)/np.sqrt(10)
print(c,RMSE)
print(c[0]+c[1]*2030)
plt.plot(t,y,'ro',t1,c[0]+c[1]*t1,'b')
plt.xlabel('Time')
plt.ylabel('daily oil production')
plt.ylabel('oil production (M bbl/day)')
plt.grid(True)
#plt.savefig('CP4o1o11c.png')
plt.show()