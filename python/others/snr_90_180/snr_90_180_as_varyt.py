import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

#### as sample

openings = np.load('snr_as90_500nm.npz')['openings']
snr90 = np.load('snr_as90_500nm.npz')['snr']
snr180 = np.load('snr_as180_500nm.npz')['snr']

f90 = interp1d(snr90,openings)

r = np.where(snr180>=snr90.min())

plt.plot(openings[r],f90(snr180[r]), 'r-', label='500 nm')


snr90 = np.load('snr_as90_5um.npz')['snr']
snr180 = np.load('snr_as180_5um.npz')['snr']

f90 = interp1d(snr90,openings)

r = np.where(snr180>=snr90.min())

plt.plot(openings[r],f90(snr180[r]),'b-',label = '5 um')


snr90 = np.load('snr_as90_50um.npz')['snr']
snr180 = np.load('snr_as180_50um.npz')['snr']

f90 = interp1d(snr90,openings)

r = np.where(snr180>=snr90.min())

plt.plot(openings[r],f90(snr180[r]),'g-',label = '50 um')


snr90 = np.load('snr_as90_500um.npz')['snr']
snr180 = np.load('snr_as180_500um.npz')['snr']

f90 = interp1d(snr90,openings)

r = np.where(snr180>=snr90.min())

plt.plot(openings[r],f90(snr180[r]),'m-',label = '500 um')

plt.xlabel(r'Semi-collection angle @ 180$^\circ$ (deg)')
plt.ylabel(r'Semi-collection angle @ 90$^\circ$  (deg)')
plt.legend(loc='upper left')
plt.show()

