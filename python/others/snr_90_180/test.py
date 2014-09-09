import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
openings = np.load('snr_wpw90_1um.npz')['openings']


snr90 = np.load('snr_as180_500nm.npz')['snr']
snr180 = np.load('snr_as180_5um.npz')['snr']


plt.plot(openings,snr90, 'r',label = 'As-Si 500 nm')
plt.plot(openings,snr180, 'b',label = 'As-Si 5 um')
plt.xlabel(r'Semi-collection angle @ 180$^\circ$ (deg)')
plt.ylabel(r'Semi-collection angle @ 90$^\circ$  (deg)')
plt.legend(loc='upper left')
plt.show()

