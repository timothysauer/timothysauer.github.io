import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../code')
import euler as euler

def ydot(t,y): return np.sin(y)
inter = [0,4]
y0 = 0. # 100.
n = 40
t,y = euler(ydot,inter,y0,n)

plt.plot(t,y0)
plt.grid(True)
plt.xlabel('t')
plt.ylabel('y')
#plt.savefig('cp6_1_9a.png')
plt.show()