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

plt.figure()
ax1 = plt.gca()
ax2 = ax1.twinx()

openings = np.load('snr_wpw0_100nm.npz')['openings']
snr = np.load('snr_wpw0_100nm.npz')['snr']
p2b = np.load('snr_wpw0_100nm.npz')['p2b']

p10, = ax1.plot(openings,snr, 'r-', label=r'0$^\circ$')
p20, = ax2.plot(openings,p2b, 'r--', label=r'0$^\circ$')

openings = np.load('snr_wpw90_100nm.npz')['openings']
snr = np.load('snr_wpw90_100nm.npz')['snr']
p2b = np.load('snr_wpw90_100nm.npz')['p2b']

p11, = ax1.plot(openings,snr, 'b-', label=r'90$^\circ$')
p21, = ax2.plot(openings,p2b, 'b--', label=r'90$^\circ$')

openings = np.load('snr_wpw180_100nm.npz')['openings']
snr = np.load('snr_wpw180_100nm.npz')['snr']
p2b = np.load('snr_wpw180_100nm.npz')['p2b']

p12, = ax1.plot(openings,snr, 'g-', label=r'180$^\circ$')
p22, = ax2.plot(openings,p2b, 'g--', label=r'180$^\circ$')



ax1.set_ylabel(r'S/N (solid)')
ax2.set_ylabel(r'P/B (dashed)')
ax2.set_yscale('log')


ax1.legend([p10,p20,p11,p21,p12,p22],[r'0$^\circ$',r'0$^\circ$',r'90$^\circ$',r'90$^\circ$',r'180$^\circ$',r'180$^\circ$'],bbox_to_anchor=(.8, 1),loc=1,ncol=3)
ax1.set_xlabel(r'Collection semi-angle $\phi$ (deg)')

ax1.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
# ax2.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))

plt.show()

