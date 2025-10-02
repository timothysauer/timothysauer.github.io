import numpy as np
import matplotlib.pyplot as plt

def dftfilter(inter, x, m, n, p):
    c, d = inter
    t = c + (d - c)* np.arange(n)/n  # Time points for data (n)
    tp = c + (d - c)* np.arange(p)/p  # Time points for interpolant (p)
    
    # Compute interpolation coefficients using FFT
    y = np.fft.fft(x)  # Compute DFT of x
    yp = np.zeros(p, dtype=complex)  # Will hold coefficients for inverse FFT
    yp[:m//2] = y[:m//2]  # Keep only first m frequencies
    yp[m//2] = np.real(y[m//2])  # Keep the cos term only

    if m < n:  # Unless at the maximum frequency
        yp[p - m//2] = yp[m//2]  # Add complex conjugate to upper tier

    yp[p - m//2 + 1:p] = y[n - m//2 + 1:n]  # More conjugates for upper tier

    xp = np.real(np.fft.ifft(yp))*(p/n)  # Invert FFT to recover data

    # Plot data and least squares approximation
    plt.plot(t, x, 'o', label='Original Data')
    plt.plot(tp, xp, label='Least Squares Approximation')
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.legend()
    plt.title('Least Squares Trigonometric Fit')
    plt.show()

    return xp

# Example usage
# Setup parameters: interval, data points, and other settings
inter = [0, 2*np.pi]  # Example interval [0, 2π]
n = 64  # Number of data points
x = np.sin(2*np.pi * np.linspace(0, 1, n)) + 0.1*np.random.normal(size=n)  # Example data with noise
m = 16  # Number of fitting frequencies
p = 128  # Number of points for the interpolant, p >= n

xp = dftfilter(inter, x, m, n, p)
