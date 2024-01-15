import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft

t, x = np.loadtxt('dump/x_dump.txt', unpack=True)
dt = t[1] - t[0]

fourier_components, freqs = fft.fft(x), fft.fftfreq(len(x), dt)

# Mask out frequencies less than zero
fourier_components = fourier_components[freqs >= 0]
freqs = freqs[freqs >= 0]
fourier_components[0] = 0 # remove the amplitude peak at zero

peak_freq_ind = np.argmax(np.abs(fourier_components))
print(f'Dominant frequency = {freqs[peak_freq_ind]} inverse femtoseconds')

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 9.5))

# Original x-signal
ax1.plot(t, x, color='blue')
ax1.set_title('Original x-signal of the atom (id = 28)')
ax1.set_xlabel('Time (fs)')
ax1.set_ylabel('Amplitude (Å)')

ax2.plot(freqs, abs(fourier_components), 'r-')
ax2.set_title('Magnitude Spectrum')
ax2.set_xlabel('Frequency (fs$^{-1}$)')
ax2.set_ylabel('Amplitude (Å)')

plt.savefig('x_pos_spectrogram.png', dpi=100)
plt.show()