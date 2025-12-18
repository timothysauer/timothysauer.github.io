import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

plt.rcParams.update({'font.size': 14})

def ydot(t,y): return [y[1], 18*y[0]**2]
def F(s):
    a, b, yb = 1., 2., 1./12
    sol = solve_ivp(ydot,[a,b],[1./3,s],method='RK45')
    return sol.y[0,-1]-yb


sstar=brentq(F,-1,0)
sol = solve_ivp(ydot,[1,2],[1./3,sstar],method='RK45',t_eval=np.linspace(1, 2, 20))
print(sstar)
plt.plot(sol.t,sol.y[0])
plt.grid(True)
plt.xlabel('t')
plt.ylabel('y')
plt.savefig('cp7o1o3a.png')
plt.show()