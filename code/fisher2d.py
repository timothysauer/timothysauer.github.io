import numpy as np
import matplotlib.pyplot as plt

def fisher2d(xl, xr, yb, yt, tb, te, M, N, tsteps):
    # Define initial function
    def f(x, y):
        return 2 + np.cos(np.pi * x) * np.cos(np.pi * y)
    
    delt = (te - tb) / tsteps
    D = 1.0
    m = M + 1
    n = N + 1
    mn = m * n
    h = (xr - xl) / M
    k = (yt - yb) / N

    x = np.linspace(xl, xr, m)
    y = np.linspace(yb, yt, n)

    # Initialize w
    w = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            w[i, j] = f(x[i], y[j])

    for tstep in range(tsteps):
        v = w.reshape(mn)
        wold = w.copy()

        for it in range(3):  # Newton iterations
            b = np.zeros(mn)
            DF1 = np.zeros((mn, mn))
            DF2 = np.zeros((mn, mn))
            # Build Jacobian and RHS
            for i in range(1, m - 1):
                for j in range(1, n - 1):
                    idx = i + j * m
                    # Discrete Laplacian terms
                    DF1[idx, i - 1 + j * m] = -D / h ** 2
                    DF1[idx, i + 1 + j * m] = -D / h ** 2
                    DF1[idx, idx] = 2 * D / h ** 2 + 2 * D / k ** 2 - 1 + 1 / delt
                    DF1[idx, i + (j - 1) * m] = -D / k ** 2
                    DF1[idx, i + j * m] = 2 * w[i, j]
                    DF1[idx, i + (j +1) * m] = -D / k ** 2

                    b[idx] = -wold[i, j] / delt

            # Boundary conditions for bottom and top rows
            for i in range(m):
                # Bottom (j=0)
                idx_bottom = i
                DF1[idx_bottom, idx_bottom] = -3
                if i + 1 < m:
                    DF1[idx_bottom, i + 1] = 4
                if i + 2 < m:
                    DF1[idx_bottom, i + 2] = -1
                # Top (j=n-1)
                idx_top = i + (n - 1) * m
                DF1[idx_top, idx_top] = 3
                if i - 1 >= 0:
                    DF1[idx_top, (i -1)] = -4
                if i - 2 >= 0:
                    DF1[idx_top, (i -2)] = 1

            # Boundary conditions for left and right columns
            for j in range(n):
                # Left (i=0)
                idx_left = j * m
                DF1[idx_left, idx_left] = -3
                if j + 1 < n:
                    DF1[idx_left, (j + 1) * m] = 4
                if j + 2 < n:
                    DF1[idx_left, (j + 2) * m] = -1
                # Right (i=m-1)
                idx_right = (m - 1) + j * m
                DF1[idx_right, idx_right] = 3
                if j - 1 >= 0:
                    DF1[idx_right, (j -1) * m] = -4
                if j - 2 >= 0:
                    DF1[idx_right, (j -2) * m] = 1

            # Combine Jacobian matrices
            DF = DF1 + DF2
            # Compute F
            F = (DF1 + DF2 / 2).dot(v) + b
            # Newton step
            delta_v = np.linalg.solve(DF, F)
            v -= delta_v
            w = v.reshape(m, n)

        # Plot current solution
        plt.figure()
        X, Y = np.meshgrid(x, y, indexing='ij')
        plt.contourf(X, Y, w)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(f"Time step {
