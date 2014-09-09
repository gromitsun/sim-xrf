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

yshift = 1e-18
theta = 90
psi = 45
sa = 5

ev = np.load('spec_as_si_500nm_psi%d_theta%d_sa%d.npz' % (psi, theta, sa))['ev_arr']/1e3

y = np.load('spec_as_si_500nm_psi%d_theta%d_sa%d.npz' % (psi, theta, sa))['total']

plt.plot(ev,y+yshift,label = '500 nm')


y = np.load('spec_as_si_5um_psi%d_theta%d_sa%d.npz' % (psi, theta, sa))['total']

plt.plot(ev,y+yshift,label = '5 um')

y = np.load('spec_as_si_50um_psi%d_theta%d_sa%d.npz' % (psi, theta, sa))['total']

plt.plot(ev,y+yshift,label = '50 um')


y = np.load('spec_as_si_500um_psi%d_theta%d_sa%d.npz' % (psi, theta, sa))['total']

plt.plot(ev,y+yshift,label = '500 um')

# ax = plt.gca()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# plt.legend(ncol=1,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.legend(loc=0,ncol=2)
plt.yscale('log')
plt.ylim(1e-14,5e-6)
plt.xlim(0,13)
plt.xlabel('E (KeV)')
plt.ylabel(r'$I(E)/I_0(E_0)$')
plt.show()