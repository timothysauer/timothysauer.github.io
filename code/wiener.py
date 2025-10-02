import numpy as np
import soundfile as sf

data, samplerate = sf.read('bach.mp3')
signal_in = np.zeros(data.shape)
signal_out = np.zeros(data.shape)
for i in range(2):
    c = data[:,i]  # Work with one channel at a time
    p = 1.6  # parameter for cutoff
    noise_std = np.std(c) * 0.5  # 50 percent noise
    n = len(c)
    r = noise_std * np.random.randn(n) # Generate pure noise
    x = c + r # Add noise to the signal
    fx = np.fft.fft(x)  # Take FFT of the noisy signal
    sfx = np.conj(fx) * fx  # Power spectrum
    sfcapprox = np.maximum(sfx - n * (p * noise_std) ** 2, 0)
    phi = sfcapprox / sfx
    phi = np.nan_to_num(phi)  # Replace NaNs in phi with zero
    xout = np.real(np.fft.ifft(phi * fx)) # Reconstruct via inverse FFT
    signal_in[:,i] = x
    signal_out[:,i] = xout
sf.write('original.mp3', signal_in, samplerate)
sf.write('reduced.mp3', signal_out, samplerate)