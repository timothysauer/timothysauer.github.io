import numpy as np
import matplotlib.pyplot as plt

# Global variables for input pulse parameters
pa, pb, pulse = None, None, None

def hh(inter, ic, n):
    global pa, pb, pulse
    inp = list(map(float, input('Pulse start, end, muamps separated by spaces, e.g. 50. 51. 7.: ').strip().split()))
    pa, pb, pulse = inp[0], inp[1], inp[2]
    
    a, b = inter
    h = (b - a) / n  # Step size
    y = np.zeros((n + 1, len(ic)))  # Initialize the solution array
    y[0, :] = ic  # Enter initial conditions
    t = np.zeros(n + 1)
    t[0] = a
    
    for i in range(n):
        t[i + 1] = t[i] + h
        y[i + 1, :] = rk4step(t[i], y[i, :], h)  # Update using RK4 step

    # Plotting the input pulse
    plt.subplot(3, 1, 1)
    plt.plot([a, pa, pa, pb, pb, b], [0, 0, pulse, pulse, 0, 0])
    plt.grid()
    plt.axis([0, 100, 0, 2 * pulse])
    plt.ylabel('Input Pulse')

    # Plotting voltage
    plt.subplot(3, 1, 2)
    plt.plot(t, y[:, 0])
    plt.grid()
    plt.axis([0, 100, -100, 100])
    plt.ylabel('Voltage (mV)')

    # Plotting gating variables
    plt.subplot(3, 1, 3)
    plt.plot(t, y[:, 1], label='m')
    plt.plot(t, y[:, 2], label='n')
    plt.plot(t, y[:, 3], label='h')
    plt.grid()
    plt.axis([0, 100, 0, 1])
    plt.ylabel('Gating Variables')
    plt.legend(loc = 'lower left')
    plt.xlabel('Time (msec)')

    plt.show()
    return y

def rk4step(t, w, h):
    """ One step of the Runge-Kutta order 4 method """
    s1 = ydot(t, w)
    #print(t,h,w,s1)
    s2 = ydot(t + h / 2., w + h * s1 / 2.)
    s3 = ydot(t + h / 2., w + h * s2 / 2.)
    s4 = ydot(t + h, w + h * s3)
    return w + h * (s1 + 2. * s2 + 2. * s3 + s4) / 6.

def ydot(t, w):
    """ Hodgkin-Huxley equations """
    global pa, pb, pulse
    c = 1.
    g1 = 120.
    g2 = 36.
    g3 = 0.3
    T = (pa + pb) / 2.
    length = pb - pa
    e0 = -65.
    e1 = 50.
    e2 = -77.
    e3 = -54.4
    # Square pulse input
    in_pulse = pulse * (1 - np.sign(np.abs(t - T) - length / 2.)) / 2.
    
    v, m, n, h = w  # Unpack the variables
    z = np.zeros(4)
    
    z[0] = (in_pulse - g1 * m**3 * h * (v - e1) - g2 * n**4 * (v - e2) - g3 * (v - e3)) / c
    v = v - e0  # Adjust against resting potential
    z[1] = (1 - m) * (2.5 - 0.1 * v) / (np.exp(2.5 - 0.1 * v) - 1.) - m * (4. * np.exp(-v / 18.))
    z[2] = (1 - n) * (0.1 - 0.01 * v) / (np.exp(1 - 0.1 * v) - 1.) - n * (0.125 * np.exp(-v / 80.))
    z[3] = (1 - h) * 0.07 * np.exp(-v / 20) - h/(np.exp(3.-0.1*v)+1)
    return z
out = hh([0.,100.],[-65.0,0.,0.3,0.6],2000)                                
