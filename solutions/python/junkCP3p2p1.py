import numpy as np
import matplotlib.pyplot as plt
from nest import nest
from newtdd import newtdd

y0 = np.array([91.69,92.06,92.56,95.03,95.11,88.94,90.25,94.31,96.33,100.04])
x0 = np.array([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])+2015.
x = np.linspace(2015.,2025.,400)
print(newtdd(x0,y0))
y = nest(newtdd(x0,y0),x,x0)
print(nest(newtdd(x0,y0),2025.,x0))
plt.plot(x,y,'b',x0,y0,'ro')
plt.grid(True)
plt.xlabel('year')
plt.ylabel('oil production (M bbl/day)')
#plt.savefig('CP_3_2_1.png')
plt.show()