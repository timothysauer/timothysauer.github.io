import numpy as np
import soundfile as sf
import playsound as ps

def wiener(soundfile, noise_level, p):
    """ Program 10.3 Wiener filter
        Input:  soundfile (.wav or .mp3)
                noise_level proportion of noise
                p frequency cutoff
        Output: original and filtered soundfiles """
    data, samplerate = sf.read(soundfile)
    signal_in = np.zeros(data.shape)
    signal_out = np.zeros(data.shape)
    for i in range(2):
        c = data[:,i]  # Work with one channel at a time
        n = len(c)
        noise_std = np.std(c)*noise_level  # 50 percent noise
        x = c + noise_std*np.random.randn(n) # Add noise to the signal
        fx = np.fft.fft(x)  # Take FFT of the noisy signal
        sfx = np.conj(fx)*fx  # Power spectrum
        sfcapprox = np.maximum(sfx - n*(p*noise_std)**2, 0)
        phi = sfcapprox/sfx
        phi = np.nan_to_num(phi)  # Replace NaNs in phi with zero
        xout = np.real(np.fft.ifft(phi*fx)) # Reconstruct via inverse FFT
        signal_in[:,i] = x
        signal_out[:,i] = xout
    sf.write('original.mp3', signal_in, samplerate)
    sf.write('filtered.mp3', signal_out, samplerate)
    ps.playsound('original.mp3')
    ps.playsound('filtered.mp3')

# Example usage:
wiener('bach.mp3',0.5,1.5)