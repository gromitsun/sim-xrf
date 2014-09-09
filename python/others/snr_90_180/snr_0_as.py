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

plt.figure('snr_sa_0_as')
ax1 = plt.gca()
ax2 = ax1.twinx()

openings = np.load('snr_as0_500nm.npz')['openings']
snr = np.load('snr_as0_500nm.npz')['snr']*np.sqrt(n)
p2b = np.load('snr_as0_500nm.npz')['p2b']

p10, = ax1.plot(openings,snr, 'r-', label=r'500 nm')
p11, = ax2.plot(openings,p2b, 'b-', label='500 nm')


snr = np.load('snr_as0_5um.npz')['snr']*np.sqrt(n)
p2b = np.load('snr_as0_5um.npz')['p2b']

p20, = ax1.plot(openings,snr, 'r--', label=r'5 um')
p21, = ax2.plot(openings,p2b*10, 'b--', label=r'5 um ($\times$ 10)')

snr = np.load('snr_as0_50um.npz')['snr']*np.sqrt(n)
p2b = np.load('snr_as0_50um.npz')['p2b']

p30, = ax1.plot(openings,snr, 'r-.', label=r'50 um')
p31, = ax2.plot(openings,p2b*100, 'b-.', label=r'50 um ($\times$ 100)')


snr = np.load('snr_as0_500um.npz')['snr']*np.sqrt(n)
p2b = np.load('snr_as0_500um.npz')['p2b']

p40, = ax1.plot(openings,snr*10, 'r.', label=r'500 um ($\times$ 10)')
ax1.set_ylabel(r'S/N (red)', color = 'r')
p41, = ax2.plot(openings,p2b*1000, 'b.', label=r'500 um ($\times 10^3$)')
plt.ylabel(r'P/B (blue)', color = 'b')
# plt.yscale('log')
# plt.ylim(0,18)

ax1.legend(loc=2,ncol=2,bbox_to_anchor=(0, .83))
ax2.legend(loc=2,ncol=2)
# ax1.legend([p10,p20,p20,p21,p30,p31,p40,p41])
ax1.set_xlabel(r'Collection semi-angle $\phi$ (deg)')

# ax1.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
# ax2.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))

for tl in ax1.get_yticklabels():
	tl.set_color('r')
for tl in ax2.get_yticklabels():
	tl.set_color('b')
plt.show()

