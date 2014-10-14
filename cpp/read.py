import numpy as np
import matplotlib.pyplot as plt
import sys

from python import plt_kwargs

try:
    fname = sys.argv[1]
except IndexError:
    fname = 'out.txt'


# Read channel values
a = np.genfromtxt(fname, invalid_raise=False, skiprows=0)

# Read detector configurations
f = open(fname)
for line in f.readlines():
    if 'ev_offset' in line:
        ev_offset = float(line.split('=')[1])
    elif 'ev_gain' in line:
        ev_gain = float(line.split('=')[1])
    elif 'n_channels' in line:
        n_channels = float(line.split('=')[1])
f.close()

# Plot spectrum
ev = ev_offset + ev_gain * np.arange(n_channels)
plt.figure()
plt.title('Spectrum read from %s' % fname)
plt.plot(ev/1e3, a.T)
plt.yscale('log')
plt.ylim(plt_kwargs['ylim'])
plt.xlim(plt_kwargs['xlim'])
plt.xlabel('E (KeV)')
plt.ylabel(r'$I(E)/I_0(E_0)$')

plt.show()