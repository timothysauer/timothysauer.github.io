import numpy as np
import matplotlib.pyplot as plt
from nest import nest
from newtdd import newtdd
y0 = [15.05,14.96,15.07,14.74,14.54,13.81,13.91,13.75,13.70,13.20]
x0 = np.array([0.,5,10,15,20,25,30,35,40,45])+1980.
y0 = y0[4:]
x0 = x0[4:]
x = np.linspace(2000.,2025.,100)
y = nest(newtdd(x0,y0),x,x0)
print(nest(newtdd(x0,y0),2002.,x0))
print(nest(newtdd(x0,y0),2012.,x0))
print(nest(newtdd(x0,y0),2022.,x0))
plt.plot(x,y,'b',x0,y0,'ro')
plt.grid(True)
plt.xlabel('year')
plt.ylabel('arctic ice extent (M sq km)')
plt.savefig('AddEx3_1_2a.png')
plt.show()