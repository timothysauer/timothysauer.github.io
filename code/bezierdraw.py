import numpy as np
import matplotlib.pyplot as plt
# Enable interactive plotting with plt.ion() and plt.show().
# Run this script in an environment that supports Matplotlib and allows for interactive plotting (like a local Python installation).
# Click on the figure window to place points for the Bezier spline. To finish, click three times to define a new segment or press Enter to exit.
def bezierdraw():
    # Set up the plot
    plt.plot([-1, 1], [0, 0], 'k', [0, 0], [-1, 1], 'k')
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.title("Click to define Bezier spline points")
    plt.grid(True)
    plt.ion()  # Interactive mode on
    plt.show()
    x, y = plt.ginput(1)[0] # Get first mouse click
    p = [x, y]
    t = np.linspace(0, 1, 20)
    #while True:
    for i in range(3):
        new_p = plt.ginput(3,timeout=0,show_clicks=1) # Get three mouse clicks
        print(new_p)
        if len(new_p) < 3:
            break  # terminate if not enough p
        
        x_new, y_new = zip(*new_p)  # Extract new control p
        p.extend(zip(x_new, y_new))  # Add new points to list
        print(p)
        plt.plot([p[-4][0], p[-3][0]], [p[-4][1], p[-3][1]],'r:',markersize=5)
        plt.plot(p[-3][0], p[-3][1], 'rs')  # control point
        plt.plot([p[-2][0], p[-1][0]], [p[-2][1], p[-1][1]],'r:',markersize=5)
        plt.plot(p[-2][0], p[-2][1], 'rs')  # control point
        plt.plot(p[0][0], p[0][1], 'bo')  # starting point
        plt.plot(p[-1][0], p[-1][1], 'bo')  # ending point
        bx = 3 * (p[-3][0] - p[-4][0]) # Bezier spline equations
        by = 3 * (p[-3][1] - p[-4][1])
        cx = 3 * (p[-2][0] - p[-3][0]) - bx
        cy = 3 * (p[-2][1] - p[-3][1]) - by
        dx = p[-1][0] - p[0][0] - bx - cx
        dy = p[-1][1] - p[0][1] - by - cy
        xp = p[0][0] + t * (bx + t * (cx + t * dx))
        yp = p[0][1] + t * (by + t * (cy + t * dy))
        plt.plot(xp, yp, color='blue')
        plt.show()
        # Promote last to first and repeat
        p = p[-1][:]
        print(p)
        #p[0] = p[-1]  # Update first point with the last one
    plt.ioff()  # Interactive mode off

# Call the function to run the program
bezierdraw()
