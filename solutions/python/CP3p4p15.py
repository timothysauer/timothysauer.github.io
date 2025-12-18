import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../code')
from splinecoeff import splinecoeff
from splineplot import splineplot
from splineplotalt import splineplotalt

y0 = [15.05,14.96,15.07,14.74,14.54,13.81,13.91,13.75,13.70,13.20]
x0 = np.array([0.,5,10,15,20,25,30,35,40,45])+1980.
#y0 = np.array([91.69,92.06,92.56,95.03,95.11,88.94,90.25,94.31,96.33,100.04])
#x0 = np.array([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])+2015.
#x0 = np.array([0.,5.,10,15,20,25,30,35,40,45])+1960.
#x0 = np.array([0.,10,20,40,50,60])+1960.
x1,y1 = splineplot(x0,y0,10)
x2,y2 = splineplotalt(x0,y0,10)
#print(x1[25],y1[25])
#print(x1[64],y1[64])
#print(x1[84],y1[84])
#x2,y2 = np.array([2002.,2012.,2022.]), np.array([14.57,13.86,13.95])
plt.plot(x1,y1,'b',x0,y0,'ro',x2,y2,'g--')
plt.grid(True)
plt.xlabel('year')
plt.ylabel('arctic ice extent (M sq km)')
plt.savefig('cp3_4_15.png')
plt.show()