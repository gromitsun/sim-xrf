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

t_arr = np.load('wpw90.npz')['t_arr']
snr90 = np.load('wpw90.npz')['snr']
snr180 = np.load('wpw180.npz')['snr']
p2b90 = np.load('wpw90.npz')['p2b']
p2b180 = np.load('wpw180.npz')['p2b']

P90 = snr90**2*(p2b90+2)/p2b90
P180 = snr180**2*(p2b180+2)/p2b180

plt.figure('main')

plt.plot(t_arr*1e4, (3./snr90)**2, 'r', label = r'Vortex EM @ 90d')#$^\circ$')
plt.plot(t_arr*1e4, (3./snr180)**2, 'b', label = r'MAIA @ 180d')#$^\circ$')

# plt.plot(t_arr*1e4, 1/P90, 'r--')
# plt.plot(t_arr*1e4, 1/P180, 'b--')
# plt.plot(t_arr*1e4, 10/P90, 'r-.')
# plt.plot(t_arr*1e4, 10/P180, 'b-.')




ax=plt.gca()
ax.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))
plt.xlabel('Water layers thickness (t/um)')
plt.ylabel(r'Minimum required incident photons (n)')# ($\bar{n}$)')
plt.legend(loc='upper left')

# # # inset

t_arr = np.load('wpw90_thin.npz')['t_arr']
snr90 = np.load('wpw90_thin.npz')['snr']
snr180 = np.load('wpw180_thin.npz')['snr']
p2b90 = np.load('wpw90_thin.npz')['p2b']
p2b180 = np.load('wpw180_thin.npz')['p2b']

P90 = snr90**2*(p2b90+2)/p2b90
P180 = snr180**2*(p2b180+2)/p2b180

plt.figure('inset',figsize=(2,2.5))

plt.plot(t_arr*1e4, (3./snr90)**2, 'r', label = r'Vortex EM @ 90$^\circ$')
plt.plot(t_arr*1e4, (3./snr180)**2, 'b', label = r'MAIA @ 180$^\circ$')

# plt.plot(t_arr*1e4, 1/P90, 'r--')
# plt.plot(t_arr*1e4, 1/P180, 'b--')
# plt.plot(t_arr*1e4, 10/P90, 'r-.')
# plt.plot(t_arr*1e4, 10/P180, 'b-.')


ax=plt.gca()
ax.ticklabel_format(style='sci', axis='y', scilimits=(-2,2))

# Shink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0+.1, box.y0+.1, box.width * 0.9, box.height*.9])

plt.xlabel('t/um')
plt.ylabel(r'n')# ($\bar{n}$)')
plt.xlim(0,4)
plt.ylim(0,2e11)

plt.yticks([0,1e11,2e11])
plt.xticks([0,2,4])

plt.show()
