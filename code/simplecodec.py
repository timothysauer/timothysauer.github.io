import numpy as np
import sounddevice as sd

def simplecodec(x, Fs, b, n):
    """ Program 11.1 Simple audio codec
        Input:  x input signal
                Fs sampling rate of input signal
                b number of quantization bits
                n window length
        Output: out output signal """
    L = 1  # amplitude range
    q = 2*L/(2**b - 1) # quantization
    len_x = len(x)
    nw = len_x // n  # number of windows
    x = x[:n*nw]  # truncate x to fit number of windows
    x = 0.3 * x / np.max(np.abs(x)) # Normalize input signal
    M = np.zeros((n, 2*n)) # Generate MDCT matrix
    for i in range(n):
        for j in range(2 * n):
            M[i, j] = np.cos((i + 0.5)*(j + 0.5 + n/2)*np.pi/n)
    M *= np.sqrt(2/n)
    N = M.T  # inverse MDCT
    out = []
    w = np.zeros((nw-1,2*n))
    for k in range(nw - 1):    
        x0 = x[k*n:(k+2)*n]   # Extract current window
        y0 = M @ x0           # Apply MDCT
        y1 = np.round(y0/q) # Quantize
        # Transmit or save here
        y2 = y1*q                 # Dequantize
        w[k,:] = N @ y2            # Inverse MDCT
        if k > 0:
            w2 = w[k-1,n:2*n]
            w3 = w[k,:n]
            out.extend((w2 + w3)/2.)
    out = np.array(out)
    sd.play(0.3*x/np.max(np.abs(x)), Fs)  # Play input
    sd.wait()
    sd.play(0.3*out/np.max(np.abs(out)), Fs) # Play output
    sd.wait()
    return out

# Example usage:
t = np.linspace(0, 1, int(2**13), endpoint=False)
simplecodec(np.cos(2*np.pi*440*t), 2**13, 3, 2**5)