import numpy as np
import matplotlib.pyplot as plt

def poissonfem(xl, xr, yb, yt, M, N):
    # Define input functions
    def f(x, y):
        return 0

    def r(x, y):
        return 0

    def g1(x):
        return np.log(x ** 2 + 1)

    def g2(x):
        return np.log(x ** 2 + 4)

    def g3(y):
        return 2 * np.log(y)

    def g4(y):
        return np.log(y ** 2 + 1)

    m = M + 1
    n = N + 1
    mn = m * n

    h = (xr - xl) / M
    k = (yt - yb) / N

    x = np.linspace(xl, xr, m)
    y = np.linspace(yb, yt, n)

    A = np.zeros((mn, mn))
    b = np.zeros(mn)

    # Assembly process
    for i in range(1, m - 1):  # interior points
        for j in range(1, n - 1):
            xi = x[i]
            yj = y[j]

            # Compute sum of r at neighboring points
            rsum = (
                r(xi - 2 * h / 3, yj - k / 3)
                + r(xi - h / 3, yj - 2 * k / 3)
                + r(xi + h / 3, yj - k / 3)
                + r(xi + 2 * h / 3, yj + k / 3)
                + r(xi + h / 3, yj + 2 * k / 3)
                + r(xi - h / 3, yj + k / 3)
            )

            row = i + j * m

            A[row, row] = 2 * (h ** 2 + k ** 2) / (h * k) - (h * k) * rsum / 18

            # Off-diagonal entries
            A[row, (i - 1) + j * m] = -k / h - (h * k) * (
                r(xi - h / 3, yj + k / 3) + r(xi - 2 * h / 3, yj - k / 3)
            ) / 18

            A[row, (i - 1) + (j - 1) * m] = - (h * k) * (
                r(xi - 2 * h / 3, yj - k / 3) + r(xi - h / 3, yj - 2 * k / 3)
            ) / 18

            A[row, i + (j - 1) * m] = -h / k - (h * k) * (
                r(xi - h / 3, yj - 2 * k / 3) + r(xi + h / 3, yj - k / 3)
            ) / 18

            A[row, (i + 1) + (j - 1) * m] = -k / h - (h * k) * (
                r(xi + h / 3, yj - k / 3) + r(xi + 2 * h / 3, yj + k / 3)
            ) / 18

            A[row, (i + 1) + j * m] = - (h * k) * (
                r(xi + 2 * h / 3, yj + k / 3) + r(xi + h / 3, yj + 2 * k / 3)
            ) / 18

            A[row, i + (j + 1) * m] = -h / k - (h * k) * (
                r(xi + h / 3, yj + 2 * k / 3) + r(xi - h / 3, yj + k / 3)
            ) / 18

            A[row, (i + 1) + (j + 1) * m] = - (h * k) * (
                r(xi + 2 * h / 3, yj + k / 3) + r(xi + h / 3, yj + 2 * k / 3)
            ) / 18

            # Compute sum of f at neighboring
