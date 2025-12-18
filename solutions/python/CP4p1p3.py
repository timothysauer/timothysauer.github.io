import numpy as np
import matplotlib.pyplot as plt

y = [3015470894.,3694683794.,4447606236.,6171702993.,7021732148.,7887001292.]
t = np.array([0.,10.,20.,40.,50.,60.])+1960.
t1 = np.linspace(1960.,2020.,100)

A = np.array([np.ones((6)),t,t**2]).T
c = np.linalg.solve(A.T@A,A.T@y)
RMSE = np.linalg.norm(y-A@c)/np.sqrt(6)
print(c,RMSE)
print(c[0]+c[1]*1990+c[2]*1990**2)
plt.plot(t,y,'ro',t1,c[0]+c[1]*t1+c[2]*t1**2,'b')
plt.xlabel('year')
plt.ylabel('population')
plt.grid(True)
plt.savefig('CP4o1o3b.png')
plt.show()