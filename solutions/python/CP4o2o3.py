import numpy as np
import matplotlib.pyplot as plt

y = [3015470894.,3694683794.,4447606236.,6171702993.,7021732148.,7887001292.]
t = np.array([0.,10.,20.,40.,50.,60.])
t1 = np.linspace(0.,60.,100)

A = np.array([np.ones((6)),t]).T
c = np.linalg.solve(A.T@A,A.T@np.log(y))
RMSE = np.linalg.norm(y-A@c)/np.sqrt(6)
print(c,RMSE)
print(c[0]+c[1]*30)
plt.plot(t,y,'ro',t1,np.exp(c[0])*np.exp(c[1]*t1),'b')
plt.xlabel('year')
plt.ylabel('population')
plt.grid(True)
plt.savefig('CP4o2o3.png')
plt.show()