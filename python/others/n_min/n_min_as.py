import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
# from matplotlib.ticker import LinearLocator, FormatStrFormatter
# from matplotlib.colors import LogNorm

RECALC = True

if RECALC:


	from scipy.integrate import dblquad
	import sys, os
	sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/../../core/')
	from sim_calcspec import *

	protein = compound(CP = 'C30H50O10N9S')
	water = compound(CP = 'H2O')
	Zn = compound(CP = 'Zn')
	Si = compound(CP = 'Si')

	# p1 = mix_compound([water, protein],[.75,.15])
	p1 = Si

	ev = 1.2e4
	Z = symbol2number('As')
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
		
	np.savez("scatt_bg3_as.npz",ray=ray_sum,comp=comp_sum,xrf=xrf,subtend=subtend_arr)

ray=np.load("scatt_bg3_as.npz")['ray']
comp=np.load("scatt_bg3_as.npz")['comp']
xrf=np.load("scatt_bg3_as.npz")['xrf']
subtend=np.load("scatt_bg3_as.npz")['subtend']

matrix_density = 1
element_density = 5.73
u_scatt = (ray+comp)*matrix_density*subtend
u_xrf = xrf/np.pi/4*subtend*element_density
c = 1e-3#np.logspace(-9,-3,20)
x = 9+18*.02*u_scatt/(c*u_xrf)


plt.figure()
plt.plot(x)
# plt.yscale('log')
plt.xlabel('Mass concentration')
plt.ylabel('x')
plt.show()