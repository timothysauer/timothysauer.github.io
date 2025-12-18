import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

plt.rcParams.update({'font.size': 14})

def ydot(t,y): return [y[1], 2*np.exp(-2*y[0])*(1-t**2)]
def F(s):
    a, b, yb = 0., 1., np.log(2)
    sol = solve_ivp(ydot,[a,b],[0.,s],method='RK45')
    return sol.y[0,-1]-yb


sstar=brentq(F,-1,1)
sol = solve_ivp(ydot,[0,1.],[0.,sstar],method='RK45',t_eval=np.linspace(0, 1, 20))
print(sstar)
plt.plot(sol.t,sol.y[0])
plt.grid(True)
plt.xlabel('t')
plt.ylabel('y')
plt.savefig('cp7o1o3b.png')
plt.show()