import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.integrate import dblquad
import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/../../core/')
from sim_calcspec import *

p1 = compound(CP='Si')


ev = 1.2e4

Z = 33
xrf = mac_pe_line_kissel_cascade(ev, Z, KA_LINE)


_dcs_rayleigh_pol = lambda beta_rad, theta_rad:p1.dmac_rayleigh_pol(ev, theta_rad, beta_rad)*np.sin(theta_rad)
_dcs_compton_pol = lambda beta_rad, theta_rad:p1.dmac_compton_pol(ev, theta_rad, beta_rad)*np.sin(theta_rad)

ray_sum = np.array([])
comp_sum = np.array([])
subtend_arr = np.array([])
for x in range(1,91):
	omega=solid_angle(angle_range = [90-x, 90+x, -x, x], n_theta = 4*x, n_beta = 4*x)
	x = np.radians(x)
	a,b,gfun,hfun = [np.pi/2-x, np.pi/2+x, lambda y: -x, lambda y: x]
	ray = dblquad(_dcs_rayleigh_pol,a,b,gfun,hfun)[0]/omega.subtend
	comp = dblquad(_dcs_compton_pol,a,b,gfun,hfun)[0]/omega.subtend
	ray_sum = np.append(ray_sum,ray)
	comp_sum = np.append(comp_sum,comp)
	subtend_arr = np.append(subtend_arr,omega.subtend)
	
	
im = np.array([r*xrf/(4*np.pi)*np.sqrt(subtend_arr/(r*xrf/(4*np.pi)+0.05*2*(ray_sum+comp_sum))) for r in np.logspace(-9,-3,10)])
print im.shape
plt.imshow(im,extent=[1,90,1e-9,1e-3],origin = 'lower',norm=LogNorm())
plt.colorbar()
plt.yscale('log')
plt.xlabel('Collection semi-angle (deg)')
plt.ylabel('Element concentration (wt%)')
plt.show()