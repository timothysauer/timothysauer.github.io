import numpy as np
import matplotlib.pyplot as plt

def dftinterp(inter, x, n, p):
    c, d = inter
    t = c + (d - c) * np.arange(n) / n  # Original data points
    tp = c + (d - c) * np.arange(p) / p  # Time points for interpolant
    # Apply DFT
    y = np.fft.fft(x)  # Compute the Fourier coefficients
    yp = np.zeros(p, dtype=complex)  # Allocate for inverse FFT coefficients
    # Move frequencies from n to p
    yp[:n // 2 + 1] = y[:n // 2 + 1]  # Lower half frequencies
    yp[p - n // 2 + 1:p] = y[n // 2 + 1:n]  # Upper half frequencies
    # Invert FFT to recover data
    xp = np.real(np.fft.ifft(yp)) * (p / n)
    # Plot original data points and interpolant
    plt.plot(t, x, 'o', label='Data Points')
    plt.plot(tp, xp, label='Interpolant')
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.title('Fourier Interpolation')
    plt.legend()
    plt.show()

    return xp

# Example usage
# Setup parameters for the interpolation
inter = [0, 2 * np.pi]  # Interval [0, 2π]
n = 16  # Number of data points
x = np.sin(np.linspace(0, 2 * np.pi, n))  # Data points (sine function)
p = 64  # Points for the interpolant, must be >= n

# Perform Fourier interpolation
xp = dftinterp(inter, x, n, p)
