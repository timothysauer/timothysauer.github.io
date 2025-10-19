import numpy as np
import matplotlib.pyplot as plt
from nest import nest
from newtdd import newtdd

def f(x): return np.exp(np.sin(x))
x0 = np.linspace(0.,4.,11)
y0 = f(x0)
x = np.linspace(0.,4.,400)
y = nest(newtdd(x0,y0),x,x0)
y1 = f(x)
print(np.max(np.abs(y-y1)))
plt.plot(x,y,'b',x0,y0,'ro')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('AddEx3_2_2.png')
plt.show()