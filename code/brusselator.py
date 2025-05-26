import numpy as np
import matplotlib.pyplot as plt

def brusselator(xl, xr, yb, yt, tb, te, M, N, tsteps):
    Dp = 1
    Dq = 8
    C = 4.5
    K = 9

    delt = (te - tb) / tsteps
    m = M + 1
    n = N + 1
    mn = m * n
    mn2 = 2 * mn
    h = (xr - xl) / M
    k = (yt - yb) / N
    x = np.linspace(xl, xr, m)
    y = np.linspace(yb, yt, n)

    p = np.zeros((m, n))
    q = np.zeros((m, n))

    for i in range(m):          # Define initial conditions
        for j in range(n):
            p[i, j] = C + 0.1
            q[i, j] = K / C + 0.2

    for tstep in range(tsteps):
        v = np.concatenate((p.flatten(), q.flatten()))
        pold = p.copy()
        qold = q.copy()

        for it in range(3):  # Newton iterations
            DF1 = np.zeros((mn2, mn2))
            DF3 = np.zeros((mn2, mn2))
            b = np.zeros(mn2)

            for i in range(1, m-1):
                for j in range(1, n-1):
                    DF1[i + (j - 1) * m, i - 1 + (j - 1) * m] = -Dp / h**2
                    DF1[i + (j - 1) * m, i + (j - 1) * m] = (Dp * (2 / h**2 + 2 / k**2) + K + 1 + 1 / delt)
                    DF1[i + (j - 1) * m, i + 1 + (j - 1) * m] = -Dp / h**2
                    DF1[i + (j - 1) * m, i + (j - 2) * m] = -Dp / k**2
                    DF1[i + (j - 1) * m, i + j * m] = -Dp / k**2
                    b[i + (j - 1) * m] = -pold[i, j] / delt - C

                    DF1[mn + i + (j - 1) * m, mn + i - 1 + (j - 1) * m] = -Dq / h**2
                    DF1[mn + i + (j - 1) * m, mn + i + (j - 1) * m] = (Dq * (2 / h**2 + 2 / k**2) + 1 / delt)
                    DF1[mn + i + (j - 1) * m, mn + i + 1 + (j - 1) * m] = -Dq / h**2
                    DF1[mn + i + (j - 1) * m, mn + i + (j - 2) * m] = -Dq / k**2
                    DF1[mn + i + (j - 1) * m, mn + i + j * m] = -Dq / k**2
                    DF1[mn + i + (j - 1) * m, i + (j - 1) * m] = -K

                    DF3[i + (j - 1) * m, i + (j - 1) * m] = -2 * p[i, j] * q[i, j]
                    DF3[i + (j - 1) * m, mn + i + (j - 1) * m] = -p[i, j] ** 2
                    DF3[mn + i + (j - 1) * m, i + (j - 1) * m] = 2 * p[i, j] * q[i, j]
                    DF3[mn + i + (j - 1) * m, mn + i + (j - 1) * m] = p[i, j] ** 2
                    b[mn + i + (j - 1) * m] = -qold[i, j] / delt
                    
#   IT LOOKS LIKE CHAT CRAPPED OUT HERE AND DID NOT FINISH!!!!
