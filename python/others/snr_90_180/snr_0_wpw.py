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

plt.figure('snr_sa_0_wpw')
ax1 = plt.gca()
ax2 = ax1.twinx()

openings = np.load('snr_wpw0_100nm.npz')['openings']
snr = np.load('snr_wpw0_100nm.npz')['snr']*np.sqrt(n)
p2b = np.load('snr_wpw0_100nm.npz')['p2b']

ax1.plot(openings,snr, 'r-', label=r'100 nm')
ax2.plot(openings,p2b, 'b-', label=r'100 nm')


snr = np.load('snr_wpw0_1um.npz')['snr']*np.sqrt(n)
p2b = np.load('snr_wpw0_1um.npz')['p2b']

ax1.plot(openings,snr, 'r--', label=r'1 um')
ax2.plot(openings,p2b*10, 'b--', label=r'1 um ($\times$ 10)')

snr = np.load('snr_wpw0_10um.npz')['snr']*np.sqrt(n)
p2b = np.load('snr_wpw0_10um.npz')['p2b']

ax1.plot(openings,snr, 'r-.', label=r'10 um')
ax2.plot(openings,p2b*100, 'b-.', label=r'10 um ($\times$ 100)')


snr = np.load('snr_wpw0_100um.npz')['snr']*np.sqrt(n)
p2b = np.load('snr_wpw0_100um.npz')['p2b']

ax1.plot(openings,snr*10, 'r.', label=r'100 um ($\times$ 10)')
ax1.set_ylabel(r'S/N (red)', color = 'r')
ax2.plot(openings,p2b*1000, 'b.', label=r'100 um ($\times 10^3$)')
plt.ylabel(r'P/B (blue)', color = 'b')
# plt.yscale('log')


ax1.legend(loc=2,ncol=2,bbox_to_anchor=(0, .83))
ax2.legend(loc=2,ncol=2)
ax1.set_xlabel(r'Collection semi-angle $\phi$ (deg)')

# ax1.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
# ax2.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))

for tl in ax1.get_yticklabels():
	tl.set_color('r')
for tl in ax2.get_yticklabels():
	tl.set_color('b')
plt.show()

