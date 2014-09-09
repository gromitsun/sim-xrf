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


openings = np.load('snr_wpw180_100nm.npz')['openings']
snr180 = np.load('snr_wpw180_100nm.npz')['snr']
snr0 = np.load('snr_wpw0_100nm.npz')['snr']

f180 = interp1d(snr180,openings)

r = np.where(snr0>=snr180.min())

x = openings[r]
y = f180(snr0[r])
plt.plot(y, x, 'r-', label='100 nm')


snr180 = np.load('snr_wpw180_1um.npz')['snr']
snr0 = np.load('snr_wpw0_1um.npz')['snr']

f180 = interp1d(snr180,openings)

r = np.where(snr0>=snr180.min())

x = openings[r]
y = f180(snr0[r])
plt.plot(y, x,'b-',label = '1 um')


snr180 = np.load('snr_wpw180_10um.npz')['snr']
snr0 = np.load('snr_wpw0_10um.npz')['snr']

f180 = interp1d(snr180,openings)

r = np.where(snr0>=snr180.min())

x = openings[r]
y = f180(snr0[r])
plt.plot(y, x,'g-',label = '10 um')


snr180 = np.load('snr_wpw180_100um.npz')['snr']
snr0 = np.load('snr_wpw0_100um.npz')['snr']

f180 = interp1d(snr180,openings)

r = np.where(snr0>=snr180.min())

x = openings[r]
y = f180(snr0[r])
plt.plot(y, x,'m-',label = '100 um')


#### as sample

openings = np.load('snr_as180_500nm.npz')['openings']
snr180 = np.load('snr_as180_500nm.npz')['snr']
snr0 = np.load('snr_as0_500nm.npz')['snr']

f180 = interp1d(snr180,openings)

r = np.where(snr0>=snr180.min())

x = openings[r]
y = f180(snr0[r])
plt.plot(y, x, 'r--', label='500 nm')


snr180 = np.load('snr_as180_5um.npz')['snr']
snr0 = np.load('snr_as0_5um.npz')['snr']

f180 = interp1d(snr180,openings)

r = np.where(snr0>=snr180.min())

x = openings[r]
y = f180(snr0[r])
plt.plot(y, x,'b--',label = '5 um')


snr180 = np.load('snr_as180_50um.npz')['snr']
snr0 = np.load('snr_as0_50um.npz')['snr']

f180 = interp1d(snr180,openings)

r = np.where(snr0>=snr180.min())

x = openings[r]
y = f180(snr0[r])
plt.plot(y, x,'g--',label = '50 um')


# snr180 = np.load('snr_as180_500um.npz')['snr']
# snr0 = np.load('snr_as0_500um.npz')['snr']

# f180 = interp1d(snr180,openings)

# r = np.where(snr0>=snr180.min())

# x = openings[r]
# y = f180(snr0[r])
# plt.plot(y, x,'m--',label = '500 um')

plt.ylabel(r'Collection semi-angle @ 0$^\circ$ (deg)')
plt.xlabel(r'Collection semi-angle @ 180$^\circ$  (deg)')
plt.legend(loc=0,ncol=2)
plt.show()

