import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

plt.figure()
ax1 = plt.gca()
ax2 = ax1.twinx()

openings = np.load('snr_as0_5um.npz')['openings']
snr = np.load('snr_as0_5um.npz')['snr']
p2b = np.load('snr_as0_5um.npz')['p2b']

ax1.plot(openings,snr, 'r-', label=r'0$^\circ$')
ax2.plot(openings,p2b, 'r--', label=r'0$^\circ$')

openings = np.load('snr_as90_5um.npz')['openings']
snr = np.load('snr_as90_5um.npz')['snr']
p2b = np.load('snr_as90_5um.npz')['p2b']

ax1.plot(openings,snr, 'b-', label=r'90$^\circ$')
ax2.plot(openings,p2b, 'b--', label=r'90$^\circ$')

openings = np.load('snr_as180_5um.npz')['openings']
snr = np.load('snr_as180_5um.npz')['snr']
p2b = np.load('snr_as180_5um.npz')['p2b']

ax1.plot(openings,snr, 'g-', label=r'180$^\circ$')
ax2.plot(openings,p2b, 'g--', label=r'180$^\circ$')



ax1.set_ylabel(r'S/N (solid)')
ax2.set_ylabel(r'P/B (dashed)')
ax2.set_yscale('log')


ax1.legend(loc=0,ncol=3)
ax1.set_xlabel(r'Collection semi-angle $\phi$ (deg)')

ax1.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
# ax2.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))

plt.show()

