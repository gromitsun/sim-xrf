import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.colors import LogNorm

RECALC = False

if RECALC:


	from scipy.integrate import dblquad
	import sys, os
	sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/../../core/')
	from sim_calcspec import *

	protein = compound(CP = 'C30H50O10N9S')
	water = compound(CP = 'H2O')
	Zn = compound(CP = 'Zn')

	p1 = water
	p1 = mix_compound([water, protein],[.75,.15])


	ev = 1e4
	theta = np.radians(np.arange(1,181))
	beta = np.radians(np.arange(-90,91))

	theta, beta = np.meshgrid(theta,beta)

	scatt = p1.dmac_rayleigh_pol(ev,theta,beta)+p1.dmac_compton_pol(ev,theta,beta)


	Z = symbol2number('Zn')
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
	np.savez("scatt_bg3_wpw.npz",ray=ray_sum,comp=comp_sum,xrf=xrf,subtend=subtend_arr)

data = np.load("scatt_bg3_wpw.npz")
ray=data['ray']
comp=data['comp']
xrf=data['xrf']
subtend=data['subtend']

c = np.logspace(-9,-3,10).reshape(-1,1)

im = c*xrf/(4*np.pi)*np.sqrt(subtend/(c*xrf/(4*np.pi)+0.02*2*(ray+comp)))

# # #fonts# # #
import matplotlib
from matplotlib import rc

matplotlib.rcParams['pdf.fonttype'] = 'truetype'
fontProperties = {'family':'serif','serif':['Arial'],
    'weight' : 'normal', 'size' : '12'}
rc('font',**fontProperties)
# # #

plt.figure()
plt.imshow(im,extent=[1,90,1e-9,1e-3],origin = 'lower',norm=LogNorm())
plt.colorbar()
plt.yscale('log')
plt.xlabel('Collection semi-angle (deg)')
plt.ylabel('Zn mass concentration')

plt.figure()
plt.imshow(9/(im**2),extent=[1,90,1e-9,1e-3],origin = 'lower',norm=LogNorm(), vmin = 5e1, vmax = 5e14)
plt.colorbar()
plt.yscale('log')
plt.xlabel('Collection semi-angle (deg)')
plt.ylabel('Zn mass concentration')

# # contour plot
# from matplotlib.ticker import LogFormatterMathtext
# cplot=plt.contour(9/(im**2),extent=[1,90,1e-9,1e-3],colors='w',norm=LogNorm()) #contour plot, log scale
# plt.clabel(cplot,inline=True,fontsize=15,fmt=LogFormatterMathtext()) #label on the contour, also log scale


plt.show()

X, Y = np.meshgrid(range(10),range(90))
plt.figure()
ax = plt.gca(projection='3d')
# surf = ax.plot_surface(X, Y, im, cmap=cm.coolwarm)
surf = ax.plot_surface(X, Y, np.log10(im[::-1].T), rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=True)

		
plt.yticks([0,45,90],[0,45,90])
plt.xticks([0,3,6,9],[r'10$^{-3}$',r'10$^{-5}$',r'10$^{-7}$',r'10$^{-9}$'])
plt.xlabel('Zn mass concentration')
plt.ylabel(r'Collection semi-angle $\phi$ (deg)')
ax.set_zlabel(r'Contrast parameter $\Theta$')
zlabels = {}
# for z in np.arange(-7,1,1):
	# zlabels[z+7] = r''
for z in np.arange(-7,1,1):
	zlabels[z+7] = r'10$^{%d}$'%z
ax.set_zticklabels(zlabels)
cbar = plt.colorbar(surf, shrink=0.8, aspect=10)
cmax,cmin = (np.log10(im.max()),np.log10(im.min()))
cbar.ax.get_yaxis().set_ticks([])
for j, lab in enumerate([r'10$^{-6}$',r'10$^{-5}$',r'10$^{-4}$',r'10$^{-3}$',r'10$^{-2}$',r'10$^{-1}$']):
    cbar.ax.text(1.5, (j-6-cmin) / (cmax-cmin), lab, ha='center', va='center')
cbar.ax.get_yaxis().labelpad = 15
ax1 = cbar.ax.twiny()
ax1.set_xlabel(r'$\Theta$', rotation=0)
ax1.set_xticks([])
plt.show()