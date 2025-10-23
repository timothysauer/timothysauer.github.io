import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../code')

t = np.linspace(0,45,10)
y = np.array([15.05,14.96,15.07,14.74,14.54,13.81,13.91,13.75,13.70,13.20])
A = np.array([np.ones((10)),t]).T
ATA = A.T@A
ATy = A.T@y
c = np.linalg.solve(ATA,ATy)
print(ATA)
print(ATy)
print(c)
print(c[0]+c[1]*50.)
plt.plot(t+1980,y,'ro',t+1980,c[0]+c[1]*t,'b')
plt.xlabel('Time')
plt.ylabel('ice extent (M sq km)')
plt.grid(True)
plt.savefig('addex4o1o2.png')
plt.show()

