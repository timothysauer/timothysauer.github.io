import numpy as np
import matplotlib.pyplot as plt

def dftinterp(inter, x, p):
    """ Program 10.1 Fourier Interpolation
        Input:  [c,d] interval
                x data points
                p length of interpolant (even number >= len(x))"""
    c, d = inter
    n = x.shape[0]
    t = c + (d - c)*np.arange(n)/n  # Original data points
    tp = c + (d - c)*np.arange(p)/p  # Time points for interpolant
    y = np.fft.fft(x)  # Compute the Fourier coefficients
    yp = np.zeros(p, dtype=complex)  # Allocate for inverse FFT coefficients
    # Move frequencies from n to p
    yp[:n//2 + 1] = y[:n//2 + 1]  # Lower half frequencies
    yp[p - n//2 + 1:p] = y[n//2 + 1:n]  # Upper half frequencies
    xp = np.real(np.fft.ifft(yp))*(p/n) # Invert FFT to recover data
    return tp, xp

# Example usage
tp, xp = dftinterp([0,2*np.pi],np.sin(np.arange(16)*2*np.pi/16),64)
plt.plot(t, x, 'o', label='Data Points') # Plot original data 
plt.plot(tp, xp, label='Interpolant') # Plot interpolant
plt.xlabel('Time')
plt.ylabel('Value')
plt.title('Fourier Interpolation')
plt.legend()
plt.show()