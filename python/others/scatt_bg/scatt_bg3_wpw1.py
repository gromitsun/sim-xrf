import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.colors import LogNorm

RECALC = True

if RECALC:


	from scipy.integrate import dblquad
	import sys, os
	sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/../../core/')
	from sim_calcspec import *

	protein = compound(CP = 'C30H50O10N9S')
	water = compound(CP = 'H2O')
	Zn = compound(CP = 'Zn')

	p1 = mix_compound([water, protein],[.75,.15])


	ev = 1e4

	Z = symbol2number('Zn')
	xrf = mac_pe_line_kissel_cascade(ev, Z, KA_LINE)
	subtend_arr = np.array([])
	for x in range(1,91):
		omega=solid_angle(angle_range = [90-x, 90+x, -x, x], n_theta = 4*x, n_beta = 4*x)
		subtend_arr = np.append(subtend_arr,omega.subtend)
		
	im = np.array([r*xrf/(4*np.pi)*subtend_arr for r in np.logspace(-9,-3,10)])

snr =np.load("scatt_bg3_wpw.npz")['im']


	
mt = 1e-20

min1 = 1./(mt*im)
min2 = 9/(np.sqrt(mt)*snr)**2
	
plt.imshow(np.maximum(min1,min2),extent=[1,90,1e-9,1e-3],origin = 'lower',norm=LogNorm())
plt.colorbar()
plt.yscale('log')
plt.xlabel('Collection semi-angle (deg)')
plt.ylabel('Zn mass concentration')
plt.show()