import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

plt.rcParams.update({'font.size': 14})

def ydot(t,y): return [y[0]-3*y[0]*y[1], -6*(t*y[1]+np.log(y[0]))]
def F(s):
    a, b, yb = 0., 1., -2/3.
    sol = solve_ivp(ydot,[a,b],[1.,s],method='RK45',rtol=1e-12)
    return sol.y[1,-1]-yb

sstar=brentq(F,0,1)
sol = solve_ivp(ydot,[0.,1.],[1.,sstar],method='RK45',rtol=1e-10,t_eval=np.linspace(0, 1, 20))
print(sstar)
plt.plot(sol.t,sol.y[0],sol.t,sol.y[1])
plt.grid(True)
plt.xlabel('t')
plt.ylabel('y')
plt.savefig('cp7o1o5b.png')
plt.show()