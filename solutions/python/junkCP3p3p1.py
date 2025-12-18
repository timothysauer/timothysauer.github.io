import numpy as np
import matplotlib.pyplot as plt
from nest import nest
from newtdd import newtdd

ch = np.cos((np.linspace(0.,3.,4)*2.+1)*np.pi/(2*4.))
b = np.pi/4.+np.pi/4.*np.array(ch) # Base points
yb = np.cos(b)  # Corresponding cos values
c = newtdd(b, yb)  # Get coefficients 
x = np.linspace(-2.,2.,100)
y0 = np.cos(x)
y = nest(c,x,b)
plt.plot(x,y,'b',b,yb,'ro',x,y0,'k')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.show()

