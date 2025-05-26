import numpy as np
import matplotlib.pyplot as plt
# Enabled interactive plotting with plt.ion() and plt.show().
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

    # Get first mouse click
    x, y = plt.ginput(1)[0]
    points = [(x, y)]
    t = np.linspace(0, 1, 100)

    while True:
        # Get three mouse clicks
        new_points = plt.ginput(3)
        if len(new_points) < 3:
            break  # terminate if not enough points
        
        # Extract the new control points
        x_new, y_new = zip(*new_points)
        points.extend(zip(x_new, y_new))  # add new points to list

        # Plot control points
        plt.plot([points[-4][0], points[-3][0]], [points[-4][1], points[-3][1]], 'r:', markersize=5)
        plt.plot(points[-3][0], points[-3][1], 'rs')  # control point
        plt.plot([points[-2][0], points[-1][0]], [points[-2][1], points[-1][1]], 'r:', markersize=5)
        plt.plot(points[-1][0], points[-1][1], 'rs')  # spline end point
        plt.plot(points[0][0], points[0][1], 'bo')  # starting point
        plt.plot(points[-1][0], points[-1][1], 'bo')  # ending point

        # Bezier spline equations
        bx = 3 * (points[-3][0] - points[-4][0])
        by = 3 * (points[-3][1] - points[-4][1])
        cx = 3 * (points[-2][0] - points[-3][0]) - bx
        cy = 3 * (points[-2][1] - points[-3][1]) - by
        dx = points[-1][0] - points[0][0] - bx - cx
        dy = points[-1][1] - points[0][1] - by - cy

        # Calculate spline points using Horner's method
        xp = points[0][0] + t * (bx + t * (cx + t * dx))
        yp = points[0][1] + t * (by + t * (cy + t * dy))

        # Plot the Bezier spline curve
        plt.plot(xp, yp)

        # Promote last to first and repeat
        points[0] = points[-1]  # Update first point with the last one

    plt.ioff()  # Interactive mode off
    plt.show()

# Call the function to run the program
bezierdraw()
