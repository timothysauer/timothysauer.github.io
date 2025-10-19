<!DOCTYPE html>
<html>
<head>
<title>Numerical Analysis 3rd Edition Sauer </title>
<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-MML-AM_CHTML"> </script>
</head>


<h1>Additional Examples 3.4 </h1>

<b>1</b> Find \(c_1\) and \(b_3\) in the cubic spline
\[ S(x) = \left\{ \begin{array}{ll}
4+12x+c_1x^2+x^3 & \mbox{ on \([0,1]\)} \\
10+(x-1)-4(x-1)^2+(x-1)^3 & \mbox{ on \([1,3]\)} \\
4+b_3(x-3)+2(x-3)^2+(x-3)^3 & \mbox{ on \([3,5]\)}
\end{array} \right.\]
Is the spline not-a-knot?

<hr>
The sections of the splines and derivatives are:
\begin{eqnarray*}
S_1(x) &=& 4+12x+c_1x^2+x^3\\
S_1'(x) &=& 12+2c_1x+3x^2\\
S_1''(x) &=& 2c_1+6x\\ \\
S_2(x) &=& 10+(x-1)-4(x-1)^2+(x-1)^3\\
S_2'(x) &=& 1-8(x-1)+3(x-1)^2\\
S_2''(x) &=& -8+6(x-1)\\ \\
S_3(x) &=& 4+b_3(x-3)+2(x-3)^2+(x-3)^3\\
S_3'(x) &=& b_3+4(x-3)+3(x-3)^2\\
S_3''(x) &=& 4+6(x-3)
\end{eqnarray*}
Equating \(S_1''(1)=S_2''(1)\) yields \(2c_1+6 = -8,\) or \(c_1= -7.\) Equating \(S_2'(3) = S_3'(3)\) yields \(1-8(2)+3(2)^2 = b_3,\) or \(b_3 = -3.\)

<p>

Since \(S_1'''(x) = S_2'''(x) = S_3'''(x) = 6,\) it is a not-a-knot cubic spline. 

<hr><hr>

<b>2</b> Use <tt>splinecoeff.py</tt> and <tt>splineplot.py</tt> to plot the natural cubic spline through the ice extent data points in Additional Example 3.1.2.
Compare the spline estimates for 2002, 2012 and 2022 with the exact values.
<hr>
Python code for the spline:
<pre>
import numpy as np
import matplotlib.pyplot as plt
from splinecoeff import splinecoeff
from splineplot import splineplot

y0 = [15.05,14.96,15.07,14.74,14.54,13.81,13.91,13.75,13.70,13.20]
x0 = np.array([0.,5,10,15,20,25,30,35,40,45])+1980.
x1,y1 = splineplot(x0,y0,10)
x2,y2 = np.array([2002.,2012.,2022.]), np.array([14.57,13.86,13.95])
print(x1[44],y1[44])
print(x1[64],y1[64])
print(x1[84],y1[84])
plt.plot(x1,y1,'b',x0,y0,'ro',x2,y2,'ks')
plt.grid(True)
plt.xlabel('year')
plt.ylabel('arctic ice extent (M sq km)')
#plt.savefig('AddEx3_4_2.png')
plt.show()
</pre>
The plot below shows the spline. Evaluated at 2002, the spline value is 14.24, compared with the actual value 14.57. Evaluated at 2012, the spline value is 13.88, compared with the actual value 13.86.
 Evaluated at 2022, the spline value is 13.55, compared with the actual value 13.95.
<p>

<img width=400 async src="figs/AddEx3_4_2.png" align=left>

</html>
