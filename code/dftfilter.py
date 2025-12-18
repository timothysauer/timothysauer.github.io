import numpy as np
import matplotlib.pyplot as plt

def dftfilter(inter, x, m, p):
    """ Program 10.2 Fourier Filtering
        Least squares fit of n data points on interval using    
           m >= 2 even number of frequencies
        Input:  inter interval [c,d]
                x data values of length an even number
                m >= 2 number of frequencies to fit
                p even number of points to plot
        Output: xp filtered points """
    c, d = inter
    n = x.shape[0]
    t = c + (d - c)* np.arange(n)/n  # Time points for data (n)
    tp = c + (d - c)* np.arange(p)/p  # Time points for interpolant (p)
    y = np.fft.fft(x)  # Compute DFT of x
    yp = np.zeros(p, dtype=complex)  # Will hold coefficients for inverse FFT
    yp[:m//2] = y[:m//2]  # Keep only first m frequencies
    yp[m//2] = np.real(y[m//2])  # Keep the cos term only
    if m < n:  # Unless at the maximum frequency
        yp[p - m//2] = yp[m//2]  # Add complex conjugate to upper tier
    yp[p - m//2 + 1:p] = y[n - m//2 + 1:n]  # More conjugates for upper tier
    xp = np.real(np.fft.ifft(yp))*(p/n)  # Invert FFT to recover data
    return tp, xp

# Example usage
x = np.sin(np.arange(64)*2*np.pi/64) + 0.1*np.random.normal(0,1,64) 
tp, xp = dftfilter([0,2*np.pi],x, 6, 128)
plt.plot(t, x, 'o', label='Original Data')
plt.plot(tp, xp, label='Least Squares Approximation')
plt.xlabel('Time')
plt.ylabel('Value')
plt.legend()
plt.title('Least Squares Trigonometric Fit')
plt.show()