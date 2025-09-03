import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
n=500
def update(frame):
    xbob = np.cos(frame*.05)
    ybob = np.sin(frame*.05)
    bob.set_xdata([xbob])
    bob.set_ydata([ybob])
    return(bob)
fig, ax = plt.subplots(figsize=(6, 6))
bob = ax.plot([np.pi/2], [0.], 'ro', markersize=5)[0] 
ax.set_xlim([-1.2, 1.2])
ax.set_ylim([-1.2, 1.2])
ani = FuncAnimation(fig=fig, func=update, frames=n, interval = 100,repeat=False)  # Store the animation object
plt.show()
