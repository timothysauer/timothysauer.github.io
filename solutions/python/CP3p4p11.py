import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../code')
from splinecoeff import splinecoeff
from splineplot import splineplot

#y0 = [15.05,14.96,15.07,14.74,14.54,13.81,13.91,13.75,13.70,13.20]
y0 = [3015470894,3694683794,4447606236,6171702993,7021732148,7887001292]
#x0 = np.array([0.,5.,10,15,20,25,30,35,40,45])+1960.
x0 = np.array([0.,10,20,40,50,60])+1960.
x1,y1 = splineplot(x0,y0,10)
print(x1[25],y1[25])
#print(x1[64],y1[64])
#print(x1[84],y1[84])
#x2,y2 = np.array([2002.,2012.,2022.]), np.array([14.57,13.86,13.95])
plt.plot(x1,y1,'b',x0,y0,'ro')
plt.grid(True)
plt.xlabel('year')
#plt.ylabel('arctic ice extent (M sq km)')
plt.ylabel('world population')
#plt.savefig('cp3_4_11.png')
plt.show()