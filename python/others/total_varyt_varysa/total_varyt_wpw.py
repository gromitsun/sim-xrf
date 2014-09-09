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

ev = np.load('spec_water_protein_water_100nm_psi%d_theta%d_sa5.npz' % (psi, theta))['ev_arr']/1e3

y = np.load('spec_water_protein_water_100nm_psi%d_theta%d_sa5.npz' % (psi, theta))['total']

plt.plot(ev,y+yshift,label = '100 nm')


y = np.load('spec_water_protein_water_1um_psi%d_theta%d_sa5.npz' % (psi, theta))['total']

plt.plot(ev,y+yshift,label = '1 um')

y = np.load('spec_water_protein_water_10um_psi%d_theta%d_sa5.npz' % (psi, theta))['total']

plt.plot(ev,y+yshift,label = '10 um')


y = np.load('spec_water_protein_water_100um_psi%d_theta%d_sa5.npz' % (psi, theta))['total']

plt.plot(ev,y+yshift,label = '100 um')

# ax = plt.gca()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# plt.legend(ncol=1,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.legend(loc=2)
plt.yscale('log')
plt.ylim(1e-15,5e-6)
plt.xlim(0,11)
plt.xlabel('E (KeV)')
plt.ylabel(r'$I(E)/I_0(E_0)$')
plt.show()