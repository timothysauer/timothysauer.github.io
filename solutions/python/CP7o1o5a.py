import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

plt.rcParams.update({'font.size': 14})

def ydot(t,y): return [-1./y[1], 2/y[0]]
def F(s):
    a, b, yb = 1., 2., 4.
    sol = solve_ivp(ydot,[a,b],[1.,s],method='RK45',rtol=1e-12)
    return sol.y[1,-1]-yb


sstar=brentq(F,1,3)
sol = solve_ivp(ydot,[1.,2.],[1.,sstar],method='RK45',rtol=1e-10,t_eval=np.linspace(1, 2, 20))
print(sstar)
plt.plot(sol.t,sol.y[0],sol.t,sol.y[1])
plt.grid(True)
plt.xlabel('t')
plt.ylabel('y')
plt.savefig('cp7o1o5a.png')
plt.show()