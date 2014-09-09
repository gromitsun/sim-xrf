import numpy as np
import matplotlib.pyplot as plt

# # #fonts# # #
import matplotlib
from matplotlib import rc

matplotlib.rcParams['pdf.fonttype'] = 'truetype'
fontProperties = {'family':'serif','serif':['Arial'],
    'weight' : 'normal', 'size' : '12'}
rc('font',**fontProperties)
# # #

t_arr = np.load('as90.npz')['t_arr']
snr90 = np.load('as90.npz')['snr']
snr180 = np.load('as180.npz')['snr']

plt.plot(t_arr*1e4, (3./snr90)**2, 'r', label = r'Vortex EM @ 90$^\circ$')
plt.plot(t_arr*1e4, (3./snr180)**2, 'b', label = r'MAIA @ 180$^\circ$')
ax=plt.gca()
ax.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
plt.xlabel('Si layer thickness (t/um)')
plt.ylabel(r'Minimum required incident photons ($\bar{n}$)')
plt.legend(loc='upper right')
plt.ylim(2e5,9e7)
plt.show()
