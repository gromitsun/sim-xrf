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
theta = 180
psi = 90
ev = np.load('spec_water_protein_water_100nm_psi%d_theta%d_sa5.npz' % (psi, theta))['ev_arr']/1e3

y = np.load('spec_water_protein_water_100nm_psi%d_theta%d_sa5.npz' % (psi, theta))['total']

plt.plot(ev,y+yshift,label = r'$\phi$ = 5$^\circ$')

y = np.load('spec_water_protein_water_100nm_psi%d_theta%d_sa15.npz' % (psi, theta))['total']

plt.plot(ev,y+yshift,label = r'$\phi$ = 15$^\circ$')

y = np.load('spec_water_protein_water_100nm_psi%d_theta%d_sa30.npz' % (psi, theta))['total']

plt.plot(ev,y+yshift,label = r'$\phi$ = 30$^\circ$')

y = np.load('spec_water_protein_water_100nm_psi%d_theta%d_sa60.npz' % (psi, theta))['total']

plt.plot(ev,y+yshift,label = r'$\phi$ = 60$^\circ$')

# ax = plt.gca()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# plt.legend(ncol=1,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.legend(loc=2,ncol=2)
plt.yscale('log')
plt.ylim(1e-15,1e-6)
plt.xlim(0,11)
plt.xlabel('E (KeV)')
plt.ylabel(r'$I(E)/I_0(E_0)$')
plt.show()