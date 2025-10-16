import numpy as np
import matplotlib as plt
import ../nest as nest
import ../newtdd as newtdd
y0 = [15.05,14.96,15.07,14.74,14.54,13.81,13.91,13.75,13.70,13.20]
x0 = [0.,5,10,15,20,25,30,35,40,45]
x = np.linspace(0.,45.,100)
y = nest(newtdd(x,y),x,x0)
plt.plot(x,y,'b',x0,y0,'ro')