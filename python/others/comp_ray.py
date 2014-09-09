from sim_cs_func import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
Z = np.arange(1,100).reshape(-1,1)
E = np.arange(1,41).reshape(1,-1)

ev = E*1e3
# theta_rad = np.radians(np.arange(0,181))
theta_rad = np.radians(180-20)



ray = dcs_rayleigh_unpol(ev, Z, theta_rad)
comp = dcs_compton_unpol(ev, Z, theta_rad)

plt.figure()
plt.title(r'$\sigma_{Compton}/\sigma_{Rayleigh}$')
plt.imshow(comp/ray,extent=[1,40,1,99],aspect = 'auto',norm=LogNorm())
plt.colorbar()
plt.xlabel('Energy (KeV)')
plt.ylabel('Z')
plt.show()