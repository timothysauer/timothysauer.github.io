import numpy as np
import matplotlib.pyplot as plt

def burgers(xl, xr, tb, te, M, N):
    alf = 5
    bet = 4
    D = 0.05
    
    # Define the functions for initial conditions and boundary conditions
    f = lambda x: 2 * D * bet * np.pi * np.sin(np.pi * x) / (alf + bet * np.cos(np.pi * x))
    l = lambda t: 0 * t  # Left boundary condition
    r = lambda t: 0 * t  # Right boundary condition

    h = (xr - xl) / M
    k = (te - tb) / N
    m = M + 1
    n = N
    sigma = D * k / (h**2)

    w = np.zeros((m, n + 1))  # Initialize solution array
    w[:, 0] = f(np.linspace(xl, xr, m))  # Set initial conditions
    w1 = w[:, 0].copy()

    for j in range(n):
        for it in range(3):  # Newton iteration
            DF1 = np.zeros((m, m))
            DF2 = np.zeros((m, m))

            DF1 = np.diag(1 + 2 * sigma * np.ones(m)) + np.diag(-sigma * np.ones(m - 1), 1) + np.diag(-sigma * np.ones(m - 1), -1)

            DF2 = np.zeros((m, m))
            DF2[1:m - 1, :] = np.diag(k * w1[2:m] / (2 * h), 1) - np.diag(k * w1[0:m - 2] / (2 * h), -1)

            DF2[1:m - 1, :] += np.diag(k * w1[1:m - 1] / (2 * h), 1) - np.diag(k * w1[1:m - 1] / (2 * h), -1)

            DF = DF1 + DF2

            F = -w[:, j] + (DF1 + DF2 / 2) @ w1  # Using Lemma 8.11

            # Dirichlet conditions for DF
            DF[0, :] = np.zeros(m)
            DF[0, 0] = 1
            DF[m-1, :] = np.zeros(m)
            DF[m-1, m-1] = 1
            
            # Dirichlet conditions for F
            F[0] = w1[0] - l(j)
            F[m-1] = w1[m-1] - r(j)

            # Solve the linear system
            w1 = w1 - np.linalg.solve(DF, F)

        # Update the solution for the next time step
        w[:, j + 1] = w1

    # Prepare data for visualization
    x = np.linspace(xl, xr, m)
    t = np.linspace(tb, te, n + 1)

    # Create a mesh plot
    X, T = np.meshgrid(x, t)
    plt.figure(figsize=(8, 6))
    plt.contourf(X, T, w.T, levels=20, cmap='viridis')
    plt.colorbar(label='w(x,t)')
    plt.xlabel('x')
    plt.ylabel('t')
    plt.title('Solution of Burgers Equation')
    plt.show()

    return w

# Example usage
w = burgers(0, 1, 0, 1, 10, 10)
