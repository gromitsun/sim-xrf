import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# # #fonts# # #
import matplotlib
from matplotlib import rc

matplotlib.rcParams['pdf.fonttype'] = 'truetype'
fontProperties = {'family':'serif','serif':['Arial'],
    'weight' : 'normal', 'size' : '12'}
rc('font',**fontProperties)
# # #

n = 1e12

plt.figure('snr_sa_90_as')
ax1 = plt.gca()
ax2 = ax1.twinx()

openings = np.load('snr_as90_500nm.npz')['openings']
snr = np.load('snr_as90_500nm.npz')['snr']*np.sqrt(n)
p2b = np.load('snr_as90_500nm.npz')['p2b']

ax1.plot(openings,snr, 'r-', label='500 nm')
ax2.plot(openings,p2b, 'b-', label='500 nm')


snr = np.load('snr_as90_5um.npz')['snr']*np.sqrt(n)
p2b = np.load('snr_as90_5um.npz')['p2b']

ax1.plot(openings,snr, 'r--', label=r'1 um')
ax2.plot(openings,p2b, 'b--', label=r'1 um')

snr = np.load('snr_as90_50um.npz')['snr']*np.sqrt(n)
p2b = np.load('snr_as90_50um.npz')['p2b']

ax1.plot(openings,snr, 'r-.', label=r'50 um')
ax2.plot(openings,p2b, 'b-.', label=r'50 um')


snr = np.load('snr_as90_500um.npz')['snr']*np.sqrt(n)
p2b = np.load('snr_as90_500um.npz')['p2b']

ax1.plot(openings,snr, 'r.', label=r'500 um')
ax1.set_ylabel(r'S/N (red)', color = 'r')
ax2.plot(openings,p2b, 'b.', label=r'500 um')
plt.ylabel(r'P/B (blue)', color = 'b')
plt.yscale('log')


ax1.legend(loc=2,ncol=2,bbox_to_anchor=(.1, 1))
ax2.legend(loc=2,ncol=2,bbox_to_anchor=(.1, .85))
ax1.set_xlabel(r'Collection semi-angle $\phi$ (deg)')

# ax1.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
# ax2.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))

for tl in ax1.get_yticklabels():
	tl.set_color('r')
for tl in ax2.get_yticklabels():
	tl.set_color('b')
plt.show()

