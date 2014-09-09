import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

plt.figure()

file = 'snr_wpw90_100nm.npz'

openings = np.load(file)['openings']
snr = np.load(file)['snr']
p2b = np.load(file)['p2b']

P = snr**2*(p2b+2)/p2b

min_p = 1./P
min_snr = 9/snr**2


plt.plot(openings, min_snr, label = r'required by S/N $\geq$ 3')
plt.plot(openings, min_p, label = r'required by signal $\geq$ 1')
plt.plot(openings, 10*min_p, label = r'required by signal $\geq$ 10')
plt.yscale('log')
plt.legend()
plt.xlabel(r'Collection semi-angle $\phi$ (deg)')
plt.ylabel(r'Minimum required incident photons ($\bar{n}$)')
plt.show()

