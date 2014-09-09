import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
openings = np.load('snr_wpw90_1um.npz')['openings']
snr90 = np.load('snr_wpw90_1um.npz')['snr']
snr180 = np.load('snr_wpw180_1um.npz')['snr']

f90 = interp1d(snr90,openings)

r = np.where(snr180>=snr90.min())

plt.plot(openings[r],f90(snr180[r]), 'r', label='protein')


snr90 = np.load('snr_as90_5um.npz')['snr']
snr180 = np.load('snr_as180_5um.npz')['snr']

f90 = interp1d(snr90,openings)

r = np.where(snr180>=snr90.min())

plt.plot(openings[r],f90(snr180[r]), 'b',label = 'As-Si')
plt.xlabel(r'Semi-collection angle @ 180$^\circ$ (deg)')
plt.ylabel(r'Semi-collection angle @ 90$^\circ$  (deg)')
plt.legend(loc='upper left')
plt.show()

