import numpy as np
import matplotlib.pyplot as plt
from nest import nest
from newtdd import newtdd

y0 = [3015470894.,3694683794.,4447606236.,6171702993.,7021732148.,7887001292.]
x0 = np.array([0.,10.,20.,40.,50.,60.])+1960.
y0 = y0[:]
x0 = x0[:]
x = np.linspace(1960.,2020.,100)
y = nest(newtdd(x0,y0),x,x0)
print(nest(newtdd(x0,y0),1990.,x0))
plt.plot(x,y,'b',x0,y0,'ro')
plt.grid(True)
plt.xlabel('year')
plt.ylabel('population')
#plt.savefig('AddEx3_1_2a.png')
plt.show()