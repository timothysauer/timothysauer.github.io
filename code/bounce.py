import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
global x,y,hx,hy,low,high
n=500
hx, hy = 0.02, 0.29  # Speed
low, high = 0.0, 1.0  # Boundaries
x, y = 0.1, 0.1  # Initial positions
def update(frame):
    global x,y,hx,hy,low,high
    if x < low or x > high:  # Check for collisions and update velocity
        hx = -hx  # Reverse horizontal direction
    if y < low or y > high:
        hy = -hy  # Reverse vertical direction
    x += hx       # Update position
    y += hy
    ball.set_xdata([x])
    ball.set_ydata([y])
    if frame == n:
        plt.close(fig)
    return(ball)
fig, ax = plt.subplots(figsize=(6,6))
ball = ax.plot([0.1],[0.1],'ro', markersize=8)[0] 
ax.set_xlim([low,high])
ax.set_ylim([low,high])
ani = FuncAnimation(fig=fig, func=update, frames=n, interval = 100) 
plt.show()

    
